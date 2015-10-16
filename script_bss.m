%% SCRIPT_BSS.m
% This file contain the script which can run the entire algorithm. It is
% where the main options/parameters are set and the main fonctions are
% lauch. Refere to the manual for more information of the Algorithm
%
% # DATA CREATION : Data are created either by reading mat file or
% generated with physical relationship. see description below
% # FIRST STEP : First Baysian Sequential Simulation on rho and Rho
% # SECOND STEP : Second Baysian Sequential Simulation on gG and K
% # FLOW : Measure the flow in the field.
% # PLOTTHEM : Graphical visualisation of data
%
%
%
% * rho     : Electrical resisitivity [\omega]
% * sigma   : Electrical conductivity = 1/rho []
% * phi     : Porosity []
% * K       : Hydraulic conductivity []
% * Z.std   : Secondary variable std error (2D-grid.nx x grid.n elts)
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 29.01.2015
% * *Version:* 3.0

%clear all; t.t_ini=tic; clc profile on; profile viewer;
%profsave(profile('info'),'profile_results')
addpath(genpath('./.')) % Add folder and sub-folder to path
clc; %clear all;
addpath kriging kernel data_gen




%% DATA CREATION
% DATA_GEN function gather all possible way to creat the data. |method|
% struct allow for choising which way. See more description in file
% |data.creation.m|

% Grid size
gen.xmax = 300; %total length in unit
gen.ymax = 20; %total hight in unit

gen.scale.x = [1:9]; % define the subdivision of each scale
gen.scale.y = [1:6 6 6 6];


% Generation Method
gen.method              = 'fromRho';   
% 'Paolo':              load paolo initial model and fit it to the created grid
% 'fromK':              genreate with FFTMA a field and log transform it with the parameter defined below 
% 'fromRho':            idem

gen.samp                = 1;      % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.covar.modele        = [4 100 10 0; 1 1 1 0];
gen.covar.c             = [0.99; 0.01];
gen.mu                  = .27;
gen.std                 = .06;
gen.Rho.method            = 'R2'; %'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'


gen.Rho.grid.nx           = 300;
gen.Rho.grid.ny           = 60;
gen.Rho.elec.spacing      = 2;
gen.Rho.elec.config_max   = 6000;
gen.Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 1;

gen.plotit              = false;      % display graphic or not
gen.saveit              = true; % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'test_1';
gen.seed                = 123456;

%[grid, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);
data_generation(gen);

clear all % not used anymore




%% FIRST STEP
% Generation of the high resolution electrical conductivity (gG) from
% scarse point data (g) and grid (G). see more information in |BSGS.m| file
parm.n_realisation  = 1;
parm.scale          = 1:numel(grid);
parm.seed           = rand();
parm.neigh          = true;
parm.cstk           = true;
parm.fitvar         = 0;
parm.unit           = 'Electrical Conductivitiy';

% Saving
parm.saveit         = true;
parm.name           = gen.name;
parm.gen            = gen;

% Ploting
parm.plot.bsgs      = 0;
parm.plot.ns        = 0;
parm.plot.sb        = 0;
parm.plot.kernel    = 0;
parm.plot.fitvar    = 0;
parm.plot.krig      = 0;

% parm.k.range.min = [50 50];
% parm.k.range.max = [500 500];

tic; [gG, t] = BSGS(sigma,Sigma,sigma_true,grid,parm); t.gG_glo=toc;
% GG.d = gG{end}.m{end}; % take the last simulation of the last scale
% GG.std = G.std; % NOT GOOD
% GG.x=G.x;
% GG.y=G.y;


%% SECOND STEP
% Generation of the high resolution hydraulic conductivity (gG) from scarce
% point data (K) and electrical conductivity field. We used the log of k.
% plotit=0;
% n_realisation=1;
% K_log=K;
% K_log.d = log10(K.d); % convert in log10
% tic; [kK_log,t.kK] = BSGS(K_log,GG,K_true,grid,plotit,n_realisation);t.KKg=toc;
% kK=kK_log; % back-convert
% for i=1:size(kK_log,1)
%     for u=1:size(kK_log{i}.m,1)
%         kK{i}.m{u}=10.^kK_log{i}.m{u};
%     end
% end
% 
% clear kK_log K_log
