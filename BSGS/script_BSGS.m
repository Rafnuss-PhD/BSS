%% SCRIPT_BSS.m
% This file contain the script which can run the entire algorithm. It is
% where the main options/parameters are set and the main fonctions are
% lauch.
%
% # DATA CREATION : Data are created either by reading mat file or
% generated with physical relationship. 
% # BSS : Baysian Sequential Simulation
% # PLOT : Graphical visualisation of data
%
%
% Variable naming convension:
% * rho     : Electrical resisitivity [\omega.m]
% * sigma   : Electrical conductivity = 1/rho [mS/m]
% * phi     : Porosity [-]
% * K       : Hydraulic conductivity [m/s]
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)


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
clear all; close all
load('result-BSS/GEN-test_2017-03-30_08-49');
[Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
parm.k.covar = gen.covar;
parm.kernel=kern;

parm.unit='';
parm.nscore = 1;
parm.cstk = true;
parm.seed = 'shuffle';
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 0;
parm.plot.krig=0;

parm.k.nb = [0 0 0 0 0; 10 10 10 10 20];

parm.n_realisation  = 0;

parm.par = 0;
parm.par_n = 4;

% use the log of hyd. cond.
Klog=K;
Klog.d=log(Klog.d);


[Res, t, kern] = BSGS(Klog,Sigma,grid_gen,parm);

%% Weight

%% work on the weight
% parm.aggr.method='AB';
% parm.aggr.a=0:.1:1;%.01:.01:.19;
% parm.aggr.b=5;
% [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
% parm.n_realisation  = parm.par_n*numel(parm.aggr.A)*4;
% figure(929);hold on;
% i_pt_temp=1:grid_gen.nxy;
% fx = @(a,b,x) (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
% for i_a=1:numel(parm.aggr.a)
%     a=parm.aggr.a(i_a);
%     for i_b=1:numel(parm.aggr.b)
%         b=parm.aggr.b(i_b);
%         plot(i_pt_temp,fx(a,b,i_pt_temp./grid_gen.nxy));
%     end
% end

% filename = 'ResCst11111';
% parm.par_n = 1;
% parm.aggr.method='cst';
% parm.n_realisation  = parm.par_n;
% parm.aggr.w=1;


% filename = 'ResRAD';
% parm.aggr.method='rad';
% parm.aggr.a=[1 10 100];
% parm.aggr.b=[10 100 500 1000];
% [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
% parm.n_realisation  = parm.par_n*numel(parm.aggr.A)*4;

% rad=Sigma.rad(randperm(length(Sigma.rad)));
% figure; hold on;
% for i_a=1:numel(parm.aggr.a)
%     a=parm.aggr.a(i_a);
%     for i_b=1:numel(parm.aggr.b)
%         b=parm.aggr.b(i_b);
%         w = 1-rad.*(a-b.*linspace(0,1,numel(rad)));
%         w(w<0) = 0;
%         w(w>1) = 1;
%         plot(w)
%     end
% end

 
disp('Setup ok, Start')
[Res,~,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);
disp('Run finished')

% Sigma_d = Sigma.d(any( bsxfun(@eq,Sigma.X(:),Res.X(:)') & bsxfun(@eq,Sigma.Y(:),Res.Y(:)' ),2));
id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));

E1 = nan(sum(id),parm.n_realisation);
E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
    
XY = kern.XY;
Sigma_d = Sigma.d(:);
Res_m_ns = Res.m_ns;
Res_m = Res.m;
dens = parm.kernel.dens(:);
parfor i_realisation=1:parm.n_realisation
   gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
   E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
   E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
end

try 
    if strcmp(parm.aggr.method,'AB')
        parfor i_ab = 1:numel(parm.aggr.A)
            idAB=i_ab:numel(parm.aggr.A):parm.n_realisation;
            OF1(i_ab) = sqrt(mean(mean(E1(:,idAB),2).^2));
            OF2(i_ab) = sqrt(mean(mean(E2(:,idAB),2).^2));
        end
    elseif strcmp(parm.aggr.method,'rad')
        parfor i_x = 1:numel(parm.aggr.x)
            idx=i_x:numel(parm.aggr.x):parm.n_realisation;
            OF1(i_x) = sqrt(mean(mean(E1(:,idx),2).^2));
            OF2(i_x) = sqrt(mean(mean(E2(:,idx),2).^2));
        end
    else
        OF1 = sqrt(mean(mean(E1,2).^2));
        OF2 = sqrt(mean(mean(E2,2).^2));
    end
catch
    save(['result-BSS/' filename],'Res','parm','kern','Nscore','E1','E2')
end

save(['result-BSS/' filename],'Res','parm','kern','Nscore','E1','E2','OF1','OF2')
disp('file written')
%save('result-BSS/ResAB_A_1.mat','parm','E1','E2')
% save('C:\Users\Raphael\Desktop\to be compress\ResRAD_AB.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')






%% Calibration
% XY=[X(:),Y(:)];
% Sigma_d=Sigma.d(:);
% parm.n_realisation  = parm.par_n*1;
% parm.notify = 0;
% parm.aggr.fx = @(a,b,x) (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
% fun = @(x) OF_fx(Klog,Sigma,grid_gen,parm,id,Gamma_t,XY,jpdf,Sigma_d, x(1), x(2) );

% filename = 'ResRAD';
% parm.aggr.method='rad';
% parm.aggr.a=[1 10 100];
% parm.aggr.b=[10 100 500 1000];
% 
% parm.aggr.a=[10 50 100 200 500 1000 20000];
% parm.aggr.b=[10 50 100 200 500 1000 5000 10000];
% 
% [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
% parm.n_realisation  = parm.par_n*numel(parm.aggr.A)*4;
% 
% rad=Sigma.rad(randperm(length(Sigma.rad)));
% figure; hold on;
% for i_a=1:numel(parm.aggr.a)
%     a=parm.aggr.a(i_a);
%     for i_b=1:numel(parm.aggr.b)
%         b=parm.aggr.b(i_b);
%         w = 1-rad.*(a-b.*linspace(0,1,numel(rad)));
%         w(w<0) = 0;
%         w(w>1) = 1;
%         plot(w)
%     end
% end
% 
%  
% disp('Setup ok, Start')
% [Res,~,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);
% disp('Run finished')
% 
% 
% E1 = nan(sum(id),parm.n_realisation);
% E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
%   
% id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
% Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
% Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));
% XY = kern.XY;
% Sigma_d = Sigma.d(:);
% Res_m_ns = Res.m_ns;
% Res_m = Res.m;
% dens = parm.kernel.dens(:);
% parfor i_realisation=1:parm.n_realisation
%    gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
%    E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
%    E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
% end
% 
% try 
%     if strcmp(parm.aggr.method,'AB')
%         parfor i_ab = 1:numel(parm.aggr.A)
%             idAB=i_ab:numel(parm.aggr.A):parm.n_realisation;
%             OF1(i_ab) = sqrt(mean(mean(E1(:,idAB),2).^2));
%             OF2(i_ab) = sqrt(mean(mean(E2(:,idAB),2).^2));
%         end
%     elseif strcmp(parm.aggr.method,'rad')
%         parfor i_x = 1:numel(parm.aggr.A)
%             idx=i_x:numel(parm.aggr.A):parm.n_realisation;
%             OF1(i_x) = sqrt(mean(mean(E1(:,idx),2).^2));
%             OF2(i_x) = sqrt(mean(mean(E2(:,idx),2).^2));
%         end
%     else
%         OF1 = sqrt(mean(mean(E1,2).^2));
%         OF2 = sqrt(mean(mean(E2,2).^2));
%     end
% catch
%     save(['result-BSS/' filename],'Res','parm','kern','Nscore','E1','E2')
% end
% 
% save(['result-BSS/' filename],'Res','parm','kern','Nscore','E1','E2','OF1','OF2')
% disp('file written')
% %save('result-BSS/ResAB_A_1.mat','parm','E1','E2')
% % save('C:\Users\Raphael\Desktop\to be compress\ResRAD_AB.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')
% 
% m_ns_true = reshape(Nscore.forward(log(K_true(:))),size(K_true,1),size(K_true,2));
% gamma_x = variogram_gridded_perso(m_ns_true);
% E1_true = gamma_x(id)-Gamma_t_id;
% E2_true  = ksdensity([Sigma_d m_ns_true(:)],XY)-dens;
% OF1_true = sqrt(mean(mean(E1,2).^2));
% OF2_true = sqrt(mean(mean(E2,2).^2));



