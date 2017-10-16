clc; % clear the command line
addpath(genpath('./.'));  % Add folder and sub-folder to path
dbstop if error  % activate debug in error

%% DATA CREATION
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

% Grid size
gen.xmax = 240; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.sx = 10;
gen.sy = 7;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'fromPhi';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar(1).model      = 'exponential';
gen.covar(1).range0     = [27 2.7];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 240;
gen.Rho.grid.ny           = 20; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 2; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = true;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'test';
gen.seed                = 'default';

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);

%% Create joint-pdf
file='GEN-Run_1_2017-05-07_14-37';
files={'GEN-Run_2_2017-05-07_21-56','GEN-Run_3_2017-05-06_18-21','GEN-Run_4_2017-05-07_17-17','GEN-Run_5_2017-05-07_22-56','GEN-Run_6_2017-05-07_21-08','GEN-Run_7_2017-05-07_21-04'};

%load(['result-BSS/' file],'K_true','Sigma');

% Correct
% grid_G=gen.Rho.grid;
% for i_files = 1: numel(files)
%     data=dlmread(['Y:\BSGS\result-BSS\data_gen\Run_' num2str(i_files+1) '\IO-file\f001_res.dat']);
%     output.res=flipud(reshape(data(:,3),grid_G.ny,grid_G.nx));
%     Rho.d_raw           = flipud(output.res);
%     f                   = griddedInterpolant({grid_G.y,grid_G.x},Rho.d_raw,'nearest','nearest');
%     Rho.d               = f({grid_gen.y,grid_gen.x});
%     Sigma.d             = 1000./Rho.d;
%     Sigma.d_raw         = 1000./Rho.d_raw;
%     save(['result-BSS/' files{i_files}],'-append','Sigma');
% end

for i_files = 1: numel(files)
    load(['result-BSGS/' files{i_files}],'K_true','Sigma');
    K_true_log_list(:,i_files) = log(K_true(:));
    Sigma_d_list(:,i_files) =Sigma.d(:);
end


% Scott's rules
% dS = 3.5*std(Sigma_d_list(:))*numel(Sigma_d_list)^(-1/3);
% dK = 3.5*std(K_true_log_list(:))*numel(K_true_log_list(:))^(-1/3);
% kern.axis_sec = min(Sigma_d_list(:)):dS:max(Sigma_d_list(:));
% kern.axis_prim = min(K_true_log_list(:)):dK:max(K_true_log_list(:));

[kern.prior, kern.axis_prim] = ksdensity(K_true_log_list(:));
[~ , axis_sec] = ksdensity(Sigma_d_list(:));

[X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
kern.XY = [X(:),Y(:)];

for i_files = 1: numel(files)
    jpdf(:,:,i_files) = ksdensity([Sigma_d_list(:,i_files) K_true_log_list(:,i_files)],kern.XY);
end

for i_files = 1: numel(files)
    figure;
    imagesc( reshape(jpdf(:,:,i_files),numel(kern.axis_prim), numel(kern.axis_sec)));
end

kern.dens = reshape(mean(jpdf,3),numel(kern.axis_prim), numel(kern.axis_sec));



save(['result-BSGS/' file],'-append','kern');

%% BSGS
% Generation of the high resolution electrical conductivity (sSigma) from
% scarse electrical  data (sigma) and large scale inverted ERt (Sigma).
clear all; close all; addpath(genpath('./.'));
load('result-BSS/GEN-Run_1_2017-05-07_14-37');

Nscore = nscore(kern, struct('nscore', 1), 0); %Prim.d, kern.axis_prim, 'pchip', 'pchip', parm.plot.ns

[~,s_id]=min(bsxfun(@minus,kern.axis_sec,Sigma.d(:)).^2,[],2);
sec_pdf = kern.dens(:,s_id);
sec.pdf = bsxfun(@times, sec_pdf, 1./sum(sec_pdf));
sec.axis = Nscore.forward(kern.axis_prim);

parm.k.covar = gen.covar;
parm.k.covar.range0 = fliplr(gen.covar.range0) ./ [grid_gen.dy grid_gen.dx];

parm.saveit = false;

parm.seed_path = 'shuffle';
parm.seed_search = 'shuffle';
parm.seed_U = 'shuffle';

parm.k.wradius = 3;
parm.k.lookup = false;
parm.k.nb = 80;
parm.n_real  = 11*2;
parm.mg=1;

parm.aggr.sum = 1;
parm.aggr.T = .5;
parm.aggr.method = 'cst';

% use the log of hyd. cond.
hd = sampling_pt(struct('x',1:grid_gen.nx,'y',1:grid_gen.ny),K_true,2,0);
hd.d = Nscore.forward(hd.d);

f0=kern.prior ./ sum(kern.prior);
nx = grid_gen.nx;
ny = grid_gen.ny;


[Res] = BSS(nx,ny,hd,f0,sec,parm);

% image
figure(1); clf;
n=min([5 parm.n_real+1]);
subplot(n,1,1); imagesc(log(K_true));colorbar;
for i_real=1:n-1
    subplot(n,1,1+i_real); hold on;
    imagesc(Res(:,:,i_real));colorbar;
    scatter(hd.x,hd.y,[],hd.d,'filled')
    axis tight;
end


% variogram
figure(2); clf; hold on;
id = grid_gen.x<parm.k.covar(1).range0(1).*parm.k.wradius*3;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
plot(grid_gen.x(id), Gamma_t(id),'--k')
for i_real=1:parm.n_real
    gamma_x = variogram_gridded_perso(Res(:,:,i_real));
    plot(grid_gen.x(id), gamma_x(id))
end
gamma_x = variogram_gridded_perso(reshape(Nscore.forward(log(K_true(:))),ny,nx));
plot(grid_gen.x(id), gamma_x(id),'linewidth',2)

% Joint pdf
figure(3); clf;
subplot(1,2,1); imagesc(kern.axis_prim, kern.axis_sec ,kern.dens)
dens=nan(numel(kern.axis_prim),numel(kern.axis_sec), parm.n_real);
for i_real=1:parm.n_real
    r = Res(:,:,i_real);
    dens(:,:,i_real) = reshape(ksdensity([Sigma.d(:) Nscore.inverse(r(:))],kern.XY),numel(kern.axis_prim),numel(kern.axis_sec));
end
subplot(1,2,2); imagesc(kern.axis_prim, kern.axis_sec ,mean(dens,3))



%% Weight 
% parm.aggr.method='cst';
% parm.aggr.T = (0:.1:1)';

parm.aggr.method='step';
parm.aggr.T = (0:.1:1)';


% parm.aggr.method='linear';
% parm.aggr.T = [ .1 .9];
    
% parm.aggr.method='sigmoid';
% parm.aggr.T = [ .06 Inf ; .06 1000; .06  100; .06  50; .06  20; .06  10];


filename='ResStep0-1';
parm.aggr.sum = 1;

parm.par_n=48;
% parpool(parm.par_n)
parm.n_real  = parm.par_n*numel(parm.aggr.T);


disp('Setup ok, Start')
n=6;
Res=nan(ny,nx,parm.n_real*n);
for i=4:n
    ii = (i-1)*parm.n_real + (1:parm.n_real);
    Res(:,:, ii ) =  BSS(nx,ny,hd,f0,sec,parm);
    disp(['Sim:' num2str(i/n)])
end
disp('Run finished')



id = grid_gen.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
Gamma_t_id = Gamma_t(id);
XY = kern.XY;
Sigma_d = Sigma.d(:);
dens = kern.dens(:)./sum(kern.dens(:));

E1 = nan(sum(id),parm.n_real*n);
E2 = nan(numel(kern.dens),parm.n_real*n);


parfor i_real=1:parm.n_real*n
   r = Res(:,:,i_real);
   % gamma_x = variogram_gridded_perso(r);
   % E1(:,i_real) = gamma_x(id)-Gamma_t_id;
   dens_r =ksdensity([Sigma_d r(:)],XY);
   E2(:,i_real) = dens_r./sum(dens_r) - dens;
end
disp('Error Computed')


try
    %OF1=nan(size(parm.aggr.T,1),1);
    OF2=nan(size(parm.aggr.T,1),1);
    parfor i_t = 1:size(parm.aggr.T,1)
        ii_t=i_t:size(parm.aggr.T,1):parm.n_real;
     %   OF1(i_t) = sqrt(mean(mean(E1(:,ii_t),2).^2));
        OF2(i_t) = sqrt(mean(mean(E2(:,ii_t),2).^2));
    end
catch
    %save(['result-BSGS/' filename],'Res','parm','kern','Nscore','E1','E2')
end

save(['result-BSS/' filename],'parm','E1','E2','OF1','OF2')
disp('file written')

%% Plot

scatter(OF1,OF2,'filled')

filenames={'ResCst0-1','ResStep0-1'};

figure(1); clf; hold on;
for i=1:numel(filenames)
    load(['result-BSS/' filenames{i} '.mat'], 'OF1', 'OF2','parm')
    scatter(OF1,OF2,'filled')
    text(OF1+.001,OF2+.0002,strread(num2str(parm.aggr.T'),'%s'))
end

T=[parm.aggr.A(2:end) T];



scatter(sqrt(mean((gamma_x_s(id)-Gamma_t_id).^2)),sqrt(mean((ksdensity([Sigma.d(:) log(K_true(:))],parm.kernel.XY)-parm.kernel.dens(:)).^2)),'filled')

xlabel('OF_X'); ylabel('OF_Z');

legend(res_name)%,'Location','northoutside','Orientation','horizontal')

figure(3);clf; hold on;
for i_t = 1:size(parm.aggr.T,1)
    ii_t=i_t:size(parm.aggr.T,1):parm.n_real;
    plot(mean(E1(:,ii_t),2))
end


