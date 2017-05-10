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

load(['result-BSS/' file],'K_true','Sigma');

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

% Scott's rules
dS = 3.5*std(Sigma.d(:))*numel(Sigma.d)^(-1/3);
dK = 3.5*std(log(K_true(:)))*numel(log(K_true))^(-1/3);

kern.axis_sec = (min(Sigma.d(:))-.2*range(Sigma.d(:))):dS:(max(Sigma.d(:))+.2*range(Sigma.d(:)));
kern.axis_prim = (min(log(K_true(:)))-.2*range(log(K_true(:)))):dK:(max(log(K_true(:)))+.2*range(log(K_true(:))));
[X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
kern.XY = [X(:),Y(:)];

for i_files = 1: numel(files)
    load(['result-BSS/' files{i_files}],'K_true','Sigma');
    jpdf(:,:,i_files) = ksdensity([Sigma.d(:) log(K_true(:))],kern.XY);
end

for i_files = 1: numel(files)
    figure;
    imagesc( reshape(jpdf(:,:,i_files),numel(kern.axis_prim), numel(kern.axis_sec)));
end

kern.dens = reshape(mean(jpdf,3),numel(kern.axis_prim), numel(kern.axis_sec));



save(['result-BSS/' file],'-append','kern');

%% BSGS
% Generation of the high resolution electrical conductivity (sSigma) from
% scarse electrical  data (sigma) and large scale inverted ERt (Sigma).
clear all; close all
load('result-BSS/GEN-test_2017-03-30_08-49');

parm.k.covar = gen.covar;

parm.unit='';
parm.nscore = 1;
parm.par = 0;
parm.n_realisation  = 0;
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
% use the log of hyd. cond.
Klog=K;
Klog.d=log(Klog.d);

% work on the weight
parm.aggr.fx = @(parm,grid,i_scale,i_pt) 1;


[Res, t, kern] = BSGS(Klog,Sigma,grid_gen,parm);

%% Plot






