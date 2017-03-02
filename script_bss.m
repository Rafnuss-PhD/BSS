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
gen.xmax = 200; %total length in unit [m]
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
gen.method              = 'fromPhi';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar(1).model      = 'spherical';
gen.covar(1).range0     = [15 4];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
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
gen.name                = 'test';
gen.seed                = 'default';

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);



%% BSGS
% Generation of the high resolution electrical conductivity (sSigma) from
% scarse electrical  data (sigma) and large scale inverted ERt (Sigma).
clear all; close all
load('result-BSS/GEN-test_2017-01-11_10-59');

parm.k.covar = gen.covar;

parm.unit='';
parm.nscore = 1;
parm.par = 0;
parm.n_realisation  = 5;
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






