%% SCRIPT_BSS.m
% This file contain the scriptu which can run the entire algorithm. It is
% where the main options/parameters are set and the main fonctions are
% lauch. Refere to the manual for more information of the Algorithm
% 
% # DATA CREATION : Data are created either by reading mat file or
% generated with physical relationship. see description below
% # FIRST STEP : First Baysian Sequential Simulation on g and G
% # SECOND STEP : Second Baysian Sequential Simulation on gG and K
% # FLOW : Measure the flow in the field.
% # PLOTTHEM : Grafical visualisation of data
% 
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch) 
% * *Date:* 29.01.2015
% * *Version:* 3.0

%clear all; t.t_ini=tic; clc
%profile on;
%profile viewer;
%profsave(profile('info'),'profile_results')
addpath(genpath('./.')) % Add folder and sub-folder to path

%% DATA CREATION

%%
% Setting up the grid
xmax=240; %total length in unit
ymax=20; %total hight in unit

scale.x=[3 6 10]; % define the subdivision of each scale
scale.y=[3 6 6 ]; % ny(i)=2^scale.y(i)+1
scale.n=numel(scale.x); % number of scale, ie nb of scale

grid=cell(scale.n,1);

for i=1:scale.n;
    grid{i}.nx=2^scale.x(i)+1;
    grid{i}.ny=2^scale.y(i)+1;
    grid{i}.nxy=grid{i}.nx*grid{i}.ny;
    
    grid{i}.dx=xmax/(grid{i}.nx-1);
    grid{i}.dy=ymax/(grid{i}.ny-1);
    
    grid{i}.x=linspace(0, xmax, grid{i}.nx); 
    grid{i}.y=linspace(0, ymax, grid{i}.ny);
    %grid{i}.x=linspace(grid{i}.dx/2, xmax-grid{i}.dx/2, grid{i}.nx); 
    %grid{i}.y=linspace(grid{i}.dy/2, ymax-grid{i}.dy/2, grid{i}.ny);

    [grid{i}.X, grid{i}.Y]=meshgrid(grid{i}.x,grid{i}.y);
    grid{i}.xy=1:grid{i}.nxy;
end
clear i xmax ymax

%%
% DATA_CREATION function gather all possible way to creat the data. 
% |method| struct allow for choising which way. See more description in
% file |data.creation.m|
method.gen=1;       % Method of generation of g_true, K_true and rho_true | 1: data , 2: generation from equation
method.samp=1;      % Method of sampling of K and g | 1: borehole, 2:random
plotit=1;           % display graphic or not
[K_true, g_true, rho_true, K, g, G] = data_creation(grid{scale.n}, method, plotit);
clear plotit method % not used anymore

%% FIRST STEP 
% Generation of the high resolution electrical conductivity (gG) from
% scarse point data (g) and grid (G). see more information in |BSGS.m| file
% edit BSGS.m
tic; [gG, t.gG_sim] = BSGS(g,G,g_true,grid); t.gG_glo=toc;
GG.d = gG(:,:,1);
GG.std = G.std;
return

%% SECOND STEP
% Generation of the high resolution hydraulic conductivity (gG) from scarce
% point data (K) and electrical conductivity field. We used the log of k.
n_sim = 1; % Simulation run number
K_log=K;
K_log.d = log10(K.d);
tic; [kK_log,t.kK] = BSGS(K_log,GG,K_true,n_sim,0);t.KKg=toc;
kK=10.^kK_log;
clear kK_log K_log



%% FLOW
% Measure flow in the field
plotit=1;
flow(grid,KK,plotit)
flow(grid,K_true,plotit)
clear plotit


%% PLOTTHEM
% Plot stuff
sub_n=min([n_sim,3])+1;

figure;
subplot(sub_n,1,1);hold on;
pcolor(grid.x,grid.y,g_true); plot(grid.x(g.x),grid.y(g.y),'or');xlabel('x[m]');ylabel('y[m]');title('g True Value');shading flat;colorbar;
for i=2:sub_n
	subplot(sub_n,1,i);hold on;
	pcolor(grid.x,grid.y,gG(:,:,i-1)); title(['g Simulated field: ', i-1]);xlabel('x[m]');ylabel('y[m]');shading flat;colorbar;
end

figure;
subplot(sub_n,1,1);hold on;
pcolor(grid.x,grid.y,K_true); plot(grid.x(K.x),grid.y(K.y),'or');xlabel('x[m]');ylabel('y[m]');title('K True Value');shading flat;colorbar;
for i=2:sub_n
	subplot(sub_n,1,i);hold on;
	pcolor(grid.x,grid.y,kK(:,:,i-1)); title(['K Simulated field: ', i-1]);xlabel('x[m]');ylabel('y[m]');shading flat;colorbar;
end

figure; hold on; 
plot([1000:1000:(1000*length(t.gG)-1) grid.nx*grid.ny-length(g.d)],t.gG)
plot([1000:1000:(1000*length(t.kK)-1) grid.nx*grid.ny-length(K.d)],t.kK)
title('Time Spend');xlabel('pt simulated [m]');ylabel('Time elpase [s]');
legend('g-simulation','K simulation')