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
% Variable naming convension:
% * rho     : Electrical resisitivity [\omega.m]
% * sigma   : Electrical conductivity = 1/rho [mS/m]
% * phi     : Porosity [-]
% * K       : Hydraulic conductivity [m/s]
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 19.10.2015

addpath(genpath('./.')) % Add folder and sub-folder to path
clc; % clear all;


%% DATA CREATION
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

% Grid size
gen.xmax = 240; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.scale.x(i)-1) \times (2^gen.scale.y(i)-1)$ 
gen.scale.x = [1:10];
gen.scale.y = [1:6 6 6 6 6];

% Generation Method.
gen.method              = 'fromRho';% 'fromRho';   
% 'Paolo':              load paolo initial model and fit it to the created grid
% 'fromK':              genreate with FFTMA a field and log transform it with the parameter defined below 
% 'fromRho':            idem

% Generation parameter
gen.samp                = 1;                     % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 3;          % number of well or number of point
gen.covar.modele        = [4 160 15 0; 1 1 1 0]; % covariance structure
gen.covar.c             = [0.99; 0.01]; 
gen.mu                  = .27; % parameter of the first field. 
gen.std                 = .06;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 300;
gen.Rho.grid.ny           = 60; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 1;

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'Very_large_1millionscells';
gen.seed                = 123456;

% Run the function
data_generation(gen);
%[grid, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);



%% BSGS
% Generation of the high resolution electrical conductivity (sSigma) from
% scarse electrical  data (sigma) and large scale inverted ERt (Sigma).

parm.n_realisation  = 2;
parm.scale          = numel(grid);

parm.fitvar         = 0;
parm.covar          = gen.covar;

parm.seed           = rand();
parm.neigh          = true;
parm.cstk           = true;
parm.nscore         = true;
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

parm.k.range.min = [min(sigma.d(:))-2 min(Sigma.d(:))-2];
parm.k.range.max = [max(sigma.d(:))+2 max(Sigma.d(:))+2];

[sSigma, t] = BSGS(sigma,Sigma,sigma_true,grid,parm);







%% Compare field
fieldname = {'random_2015-10-26_14-08', 'random_2015-10-26_14-43'};
figure; hold on; axis equal
for i= 1:numel(fieldname)
    load(fieldname{i})
    data=[];
    for i_realisation =1:parm.n_realisation
        data= [data, Y{end}.m{i_realisation}(:)];
    end
    D = pdist(data');
    D_y = cmdscale(D);
    scatter(D_y(:,1),D_y(:,2))
end
legend('variable path','constant path')
