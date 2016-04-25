
%% Introduction

clc; % clear the command line
addpath(genpath('./.'));  % Add folder and sub-folder to path
dbstop if error  % activate debug in error


%% DATA CREATION
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

% Grid size
gen.xmax = 100; %total length in unit [m]
gen.ymax = 100; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.sx = 10;
gen.sy = 10;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal-Random'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'Normal-Random';

% Generation parameter
gen.samp                = 2;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 500;          % number of well or number of point
gen.covar.modele        = [4 20 20 0]; %; 1 1 1 0]; % covariance structure
gen.covar.c             = [1];% 0.01]; 
gen.mu                  = 0; % parameter of the first field. 
gen.std                 = 1;
%gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = '100x100-20x20';
gen.seed                = 123456;

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);


%% SGSIM
% Generation of the high resolution electrical conductivity (sSigma) from
% scarse electrical  data (sigma) and large scale inverted ERt (Sigma).
load('data_gen/data/GEN-100x100-20x20_2016-04-01_16-30');
parm.gen=gen;

parm.n_realisation  = 1;
parm.par = 0;

%parm.name           = 'test_grad_def';
%
% parm.k.nb_neigh = [0 0 0 0 0; 4 4 4 4 4];
%parm.plot.krig =1;
% 
% parm.fitvar         = 0;
% 
% parm.seed           = rand();
% parm.neigh          = true;
% parm.cstk           = false;
parm.nscore         = false;
% parm.unit           = 'Electrical Conductivitiy';
% 
% % Saving

% parm.name           = gen.name;
% 
% % Ploting
% parm.plot.bsgs      = 0;
% parm.plot.ns        = 0;
% parm.plot.sb        = 0;
% parm.plot.kernel    = 0;
% parm.plot.fitvar    = 0;
% parm.plot.krig        = 1;
% 
% parm.k.range.min = [min(sigma.d(:))-2 min(Sigma.d(:))-2];
% parm.k.range.max = [max(sigma.d(:))+2 max(Sigma.d(:))+2];

sigma.d=[];
sigma.x=[];
sigma.y=[];
sigma.n=0;

parm.cstk = false;
parm.saveit = false;
parm.scale = [7;7];

parm.n_realisation  = 50;
[~, t_nocstk_50, ~, ~, ~] = SGSIM(sigma,sigma_true,grid_gen,parm);


%% Time Computation
parm.par = 0;
sigma.d=[];
sigma.x=[];
sigma.y=[];
sigma.n=0;
parm.saveit = false;
parm.nscore         = false;
parm.scale = [7;7];

parm_cstk= [0 1];
parm_n_realisation=[1 5 10 50];

for parm_cstk_i=1:length(parm_cstk)
    for parm_n_realisation_i=1:length(parm_n_realisation)
        parm.cstk = parm_cstk(parm_cstk_i);
        parm.n_realisation  = parm_n_realisation(parm_n_realisation_i);
        [~, t{parm_cstk_i}{parm_n_realisation_i}, ~, ~, ~] = SGSIM(sigma,sigma_true,grid_gen,parm);
    end
end


t_fixe = 1;
t_init = 4;
t_krig = 20;
real_lim=50;
for real_i = 1:real_lim
    Y(real_i,:)=[t_fixe, t_init*real_i, t_krig, t_krig*(real_i-1)];
end
figure; hold on;
area(Y)
plot(1, t_cstk_1.global,'x','linewidth',5,'markersize',10)
plot(1, t_cstk_1.global-t_cstk_1.cstk,'x','linewidth',5,'markersize',10)
plot(1, t_nocstk_1.global,'x','linewidth',5,'markersize',10)

plot(5, t_cstk_5.global,'x','linewidth',5,'markersize',10)
plot(5, t_cstk_5.global-t_cstk_5.cstk,'x','linewidth',5,'markersize',10)
plot(5, t_nocstk_5.global,'x','linewidth',5,'markersize',10)

plot(10, t_cstk_10.global,'x','linewidth',5,'markersize',10)
plot(10, t_cstk_10.global-t_cstk_10.cstk,'x','linewidth',5,'markersize',10)
plot(10, t_nocstk_10.global,'x','linewidth',5,'markersize',10)

plot(50, t_cstk_50.global,'x','linewidth',5,'markersize',10)
plot(50, t_cstk_50.global-t_cstk_50.cstk,'x','linewidth',5,'markersize',10)
plot(50, t_nocstk_50.global,'x','linewidth',5,'markersize',10)

xlabel('Realisation')
ylabel('Time [sec]')
legend('Constant general','General','1 realisation kriging', 'kriging','Measured Data')


%% Unconditional Biais ?
load('data_gen/data/GEN-100x100-20x20_2016-04-01_16-30');
parm.gen=gen;
parm.par = 0;
parm.scale = [7;7];
parm.nscore = false;

sigma.d=[];
sigma.x=[];
sigma.y=[];
sigma.n=0;

parm.n_realisation  = 400;

parm.name = 'cond-400-noMG-noCstk';
parm.cstk = true;

SGSIM(sigma,sigma_true,grid_gen,parm);


load('result-SGSIM/SIM-uncond-400-noMG-Cstk_2016-04-04_14-18-22')
Res1 = cat(3, Res{1}.m{:});

load('result-SGSIM/SIM-uncond-400-noMG-noCstk_2016-04-04_13-48-08')
Res2 = cat(3, Res{1}.m{:});

figure(1); title('Spatial Average: Cstk (left), noCstk (right)')
subplot(1,2,1); imagesc(mean(Res1,3)); ac=caxis;
subplot(1,2,2); imagesc(mean(Res2,3)); caxis(ac)

figure(2); title('Spatial Average: Cstk (line 1), noCstk (line 2)')
subplot(2,4,1); imagesc( Res1(:,:,1)); ac=caxis;
subplot(2,4,2); imagesc( Res1(:,:,2)); caxis(ac)
subplot(2,4,3); imagesc( Res1(:,:,3)); caxis(ac)
subplot(2,4,4); imagesc( Res1(:,:,4)); caxis(ac)
subplot(2,4,5); imagesc( Res2(:,:,1)); caxis(ac)
subplot(2,4,6); imagesc( Res2(:,:,2)); caxis(ac)
subplot(2,4,7); imagesc( Res2(:,:,3)); caxis(ac)
subplot(2,4,8); imagesc( Res2(:,:,4)); caxis(ac)

figure(3); title('Spatial Extremum: Cstk (left), noCstk (right)')
subplot(1,2,1); imagesc(max(abs(Res1),[],3)); ac=caxis;
subplot(1,2,2); imagesc(max(abs(Res2),[],3)); caxis(ac)

figure(4); title('Spatial Variance: Cstk (left), noCstk (right)')
subplot(1,2,1); imagesc(var(Res1,0,3)); ac=caxis;
subplot(1,2,2); imagesc(var(Res2,0,3)); caxis(ac)

figure(5);clf; title('Histogram')
ksdensity2 = @(x) ksdensity(x(:));
subplot(1,2,1); hold on;
for i=1:size(Res1,3)
    ksdensity2(Res1(:,:,i))
end
plot(-5:.1:5,normpdf(-5:.1:5,0,1),'--k')
subplot(1,2,2); hold on;
for i=1:size(Res2,3)
    ksdensity2(Res2(:,:,i))
end
plot(-5:.1:5,normpdf(-5:.1:5,0,1),'--k','linewidth',2)

figure(6);clf; title('Histogram'); hold on;
ksdensity(Res1(:))
ksdensity(Res2(:))
plot(-5:.1:5,normpdf(-5:.1:5,0,1),'--k','linewidth',2)
%legend('Cstk', 'noCstk','Theorical')




%% Most simple exemple
% | 1 | 2 | 3 | 4 |
sigma_true=[];
sigma.d=[];
sigma.x=[];
sigma.y=[];
sigma.n=0;

grid_gen.sx = 1;
grid_gen.sy = 1;
grid_gen.x = 1:3;
grid_gen.y = 1:3;

parm=struct();
parm.gen.covar.modele = [4 1 1 0];
parm.gen.covar.c = 1;
parm.k.wradius = Inf; % force to take all point
parm.k.nb_neigh = [0 0 0 0 0; Inf Inf Inf Inf Inf];
parm.saveit = false;

parm.neigh = 0;
parm.nscore = 0;
parm.par = 0;


parm.n_realisation  = 100;

parm.cstk = true;
parm.seed = 1;
figure(2)
for i=1:100
    [Res, ~, ~, parm, ~] =SGSIM(sigma,sigma_true,grid_gen,parm);
    Res_t = reshape([Res{:}.m{:}],3,3,parm.n_realisation);
    
end
Res{end}.sim.xy_r


parm.cstk = false;
[Res2, ~, ~, parm2, ~] =SGSIM(sigma,sigma_true,grid_gen,parm);

Res_t = reshape([Res{:}.m{:}],3,3,parm.n_realisation);

Res_m = mean(Res_t,3);
Res_std = std(Res_t,[],3);
figure(1); hold on;
plot(Res_std(Res{end}.sim.xy_r))
plot(Res_m(Res{end}.sim.xy_r))



%% Conditional Biais ?
load('data_gen/data/GEN-100x100-20x20_2016-04-01_16-30');
parm.gen=gen;
parm.par = 0;
parm.scale = [7;7];
parm.nscore = false;

k=grid_gen.nsy*.05;
sigma = sampling_pt(grid_gen,sigma_true,2,k);

parm.n_realisation  = 400;

parm.name = 'cond-400-noMG-noCstk';
parm.cstk = true;

SGSIM(sigma,sigma_true,grid_gen,parm);


%% exemple of 
