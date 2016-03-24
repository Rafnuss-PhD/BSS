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
gen.sx = 8;
gen.sy = 5;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'fromRho';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar.modele        = [4 27 2.7 0; 1 1 1 0]; % covariance structure
gen.covar.c             = [0.99; 0.01]; 
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 200;
gen.Rho.grid.ny           = 15; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'SimilarToPaolo-new';
gen.seed                = 123456;

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);



%% BSGS
% Generation of the high resolution electrical conductivity (sSigma) from
% scarse electrical  data (sigma) and large scale inverted ERt (Sigma).
load('data_gen/data/GEN-SimilarToPaolo-new_2016-03-04_11-25');
parm.gen=gen;

parm.n_realisation  = 1;
parm.par = 0;

% parm.p_w = [0.5 0];   
parm.w_X = @(scale,depth) 1;
parm.w_Z = @(scale,depth) 1;
parm.name           = 'test_grad_def';
parm.scale = [6;5]
parm.k.nb_neigh = [0 0 0 0 0; 4 4 4 4 4];
parm.k.model    = [ 4 40 4 0; 1 1 1 0]; 
%parm.plot.krig =1;
% 
% parm.fitvar         = 0;
% 
% parm.seed           = rand();
% parm.neigh          = true;
% parm.cstk           = false;
% parm.nscore         = true;
% parm.unit           = 'Electrical Conductivitiy';
% 
% % Saving
% parm.saveit         = true;
% parm.name           = gen.name;
% 
% % Ploting
% parm.plot.bsgs      = 0;
% parm.plot.ns        = 0;
% parm.plot.sb        = 0;
% parm.plot.kernel    = 0;
% parm.plot.fitvar    = 0;
parm.plot.krig      = 1;
% 
% parm.k.range.min = [min(sigma.d(:))-2 min(Sigma.d(:))-2];
% parm.k.range.max = [max(sigma.d(:))+2 max(Sigma.d(:))+2];

[Y, t, kernel, k,filename] = BSGS_par(sigma,Sigma,sigma_true,grid_gen,parm);




%% Plot






