
%% Create a standard case

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
gen.method              = 'fromRho';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 3;          % number of well or number of point
gen.covar.modele        = [4 70 7 0];%; 1 1 1 0]; % covariance structure
gen.covar.c             = [0.99];%; 0.01]; 
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 256;
gen.Rho.grid.ny           = 32; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'Calibration-';
gen.seed                = 123456;

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);


%% Create simulation

clear all;clc; dbstop if error
familyname='4-Calibration/';

load('result/GEN-Calibration-_2016-03-15_13-44');
parm.gen=gen;

Sigma.d = Sigma.d*10;

parm.w_krig = @(x) .4;
parm.w_sec = @(x) .4;

parm.scale = [1:7;1:6 6];


[~,~,~,~,~, filename] = BSGS(sigma,Sigma,sigma_true,grid_gen,parm);
% filename = 'SIM-_2016-03-15_16-37-50';
OF(filename)

for i_realisation=1:parm.n_realisation
    of{i}.mean(i_realisation) = mean(Y{end}.m_ns{i_realisation}(:));
    of{i}.std(i_realisation) = std(Y{end}.m_ns{i_realisation}(:));
end
    
%% old stuff
parm.n_realisation = 20;
parm.familyname = familyname;
parm.par = 1;
scale_max = 7;
funperso = @(a,scale, scale_max) -erf( bsxfun(@times, log(scale./scale_max) , a.^(-1)' ));

figure;hold on
a=[10 4 2  1 0.5 .1];
for i=1:numel(a)
    plot(funperso(a(i),0:scale_max-1, scale_max-1),'x-')
    funperso(a(i),0:scale_max-1, scale_max-1)
end

a=1:10;
for i=1:numel(a)
    parm0=parm;
    parm0.name=['calibration_real1_a' sprintf('%2.2d',i)];
    %p_w_start = funperso(a(i),0:scale_max-1, scale_max-1);
    p_w_start = a<a(i);
    parm0.p_w = [p_w_start zeros(1,grid_gen.sx-numel(p_w_start))];
    BSGS(sigma,Sigma,sigma_true,grid_gen,parm0);
end


%% Compute the objectif funtion
clc;clear all;
folder='result-time-log-calibration/4-Calibration/ert-real20/';
list = dir([folder 'SIM*']);
n_sim=size(list,1);

for i=1:n_sim
    
    load([folder list(i).name],'parm','Y','grid')
    
    % of 1 + 2 : histogram
    for i_realisation=1:parm.n_realisation
        of{i}.mean(i_realisation) = mean(Y{end}.m_ns{i_realisation}(:));
        of{i}.std(i_realisation) = std(Y{end}.m_ns{i_realisation}(:));
    end
    
    % of3 : variogram
    for i_realisation=1:parm.n_realisation
        [~, gamma_y_y] = variogram_gridded_perso(Y{end}.m{i_realisation});
        of{i}.var{i_realisation} = gamma_y_y;
    end
    
    % of3 : ert
    dmin = parm.gen.Rho.dmin;
    grid_Rho = parm.gen.Rho.grid;
    elec = parm.gen.Rho.elec;
    
    dmin.filepath       = 'result/4-Calibration/IO-file/';
    delete([dmin.filepath '*']);
    dmin.header         = 'Forward';  % title of up to 80 characters
    dmin.job_type       = 0;
    dmin.readonly       = 0;
    
    % Interpolation

    
    for i_realisation=1:parm.n_realisation
        rho_true =  1000./Y{end}.m{i_realisation}; %sigma_true;%
        dmin.rho_min        = min(rho_true(:));
        dmin.rho_avg        = mean(rho_true(:));
        dmin.rho_max        = max(rho_true(:));
        
        f                   = griddedInterpolant({grid{end}.y,grid{end}.x},rho_true,'nearest','nearest');
        dmin.rho_true       = f({grid_Rho.y,grid_Rho.x});
        dmin.filename       = 'gtrue.dat';
        dmin.num_regions    = 1+numel(dmin.rho_true);
        
        dout                = Matlat2R2(grid_Rho,dmin,elec); % write file and run forward modeling
        of{i}.ert{i_realisation} = dout.output.pseudo;
        of{i}.ert_interp{i_realisation} = dout.output.pseudo_interp;
    end
    
    % of4 = mean( (mean([pseudo{:}]')'-parm.gen.Rho.f.output.pseudo).^2 );
    of{i}.p_w = parm.p_w;
end

% Reference variogram
load([folder list(i).name])
load('data_gen/data/GEN-SimilarToPaolo-27_2015-11-20_12-54');
[~, gamma_y_s] = variogram_gridded_perso(reshape(sigma.d,numel(unique(sigma.y)),numel(unique(sigma.x))));
var_ref=gamma_y_s;
myfun = @(x,h) semivariogram1D(h,1,x,parm.k.model(1),0);
id= grid_gen.y<parm.k.range(2)*1.2;
% Reference ert
ert_ref=parm.gen.Rho.f.output.pseudo;
ert_interp_ref=parm.gen.Rho.f.output.pseudo_interp;

save([folder 'of'],'of','var_ref','ert_ref','ert_interp_ref');


%% plot
for i=1:n_sim
    of_mean(:,i)=of{i}.mean;
    of_std(:,i)=of{i}.std;
end

for i=1:n_sim
    of_var(:,i)=mean(reshape([of{i}.var{:}],grid_gen.ny,parm.n_realisation),2)';
end

for i=1:n_sim
    of_ert(:,i)=mean(([of{i}.ert{:}]-repmat(ert_ref,1,parm.n_realisation)).^2);
end

figure; 
subplot(2,2,1); hold on;
plot([0 n_sim],[nanmean(Nscore.forward(Sigma.d(:))) nanmean(Nscore.forward(Sigma.d(:)))],'--k')
plot([0 n_sim],[nanmean(Nscore.forward(sigma_true(:))) nanmean(Nscore.forward(sigma_true(:)))],'--k')
plot([0 n_sim],[0 0],'-r')
boxplot(of_mean); ylabel('Average of normal Transform field'); axis tight
subplot(2,2,2); hold on;
plot([0 n_sim],[nanstd(Nscore.forward(Sigma.d(:))) nanstd(Nscore.forward(Sigma.d(:)))],'--k')
plot([0 n_sim],[nanstd(Nscore.forward(sigma_true(:))) nanstd(Nscore.forward(sigma_true(:)))],'--k')
plot([0 n_sim],[1 1],'-r')
boxplot(of_std); ylabel('Std dev. of normal Transform field'); axis tight
subplot(2,2,3); hold on; ylabel('vertical variogram')
plot(grid_gen.y(id), of_var(id,:))
plot(grid_gen.y(id), var_ref(id),'Linewidth',3)
plot(grid_gen.y(id),myfun(parm.k.range(2),grid_gen.y(id)),'--k','linewidth',3);
subplot(2,2,4);boxplot(of_ert); ylabel('RMSE of ERT measure')