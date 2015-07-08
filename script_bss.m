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
scale.y=[3 6 6];
scale.n=numel(scale.x); % number of scale, ie nb of scale

grid=cell(scale.n,1);

% construction of the grid at each scale
for i=1:scale.n;
    grid{i}.nx=2^scale.x(i)+1;
    grid{i}.ny=2^scale.y(i)+1;
    grid{i}.nxy=grid{i}.nx*grid{i}.ny; % total number of cells
    
    grid{i}.dx=xmax/(grid{i}.nx-1);
    grid{i}.dy=ymax/(grid{i}.ny-1);
    
    grid{i}.x=linspace(0, xmax, grid{i}.nx); % coordinate of cells center
    grid{i}.y=linspace(0, ymax, grid{i}.ny);
    grid{i}.xy=1:grid{i}.nxy;

    [grid{i}.X, grid{i}.Y]=meshgrid(grid{i}.x,grid{i}.y); % matrix coordinate
end
clear i xmax ymax

%%
% DATA_CREATION function gather all possible way to creat the data. 
% |method| struct allow for choising which way. See more description in
% file |data.creation.m|
method.gen =2;      % Method of generation of g_true, K_true and rho_true | 1: data (Paolo), 2: generation from equation
method.samp=1;      % Method of sampling of K and g | 1: borehole, 2:random
plotit=1;           % display graphic or not | 1: yes, 0: no
[K_true, g_true, rho_true, K, g, G] = data_creation(grid{scale.n}, method, plotit);
clear plotit method % not used anymore

%% FIRST STEP 
% Generation of the high resolution electrical conductivity (gG) from
% scarse point data (g) and grid (G). see more information in |BSGS.m| file
% edit BSGS.m
plotit=0;
n_sim=2;
tic; [gG, t.gG_sim] = BSGS(g,G,g_true,grid,plotit,n_sim); t.gG_glo=toc;
GG.d = gG{end}.m;
GG.std = G.std; % NOT GOOD
GG.x=G.x;
GG.y=G.y;
return

%% SECOND STEP
% Generation of the high resolution hydraulic conductivity (gG) from scarce
% point data (K) and electrical conductivity field. We used the log of k.
K_log=K;
K_log.d = log10(K.d);
tic; [kK_log,t.kK] = BSGS(K_log,GG,K_true,grid);t.KKg=toc;
kK=kK_log;
for i=1:size(kK_log,1)
    kK{i}.m=10.^kK_log{i}.m;
end

clear kK_log K_log



%% PLOTTHEM
% Plot stuff
sub_n=min([n_sim,5])+1;

if sub_n<4, k=sub_n;m=1;
else k=ceil(sub_n/2); m=2; end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(k,m,1); hold on; xlabel('x[m]');ylabel('y[m]');
pcolor(grid{end}.x,grid{end}.y,g_true);shading flat;colorbar; caxis([min(g_true(:)) max(g_true(:))]); axis tight
title([ 'True Electrical Conductivity g_{true}: \mu=' num2str(mean2(g_true)) ' | \sigma='   num2str(std2(g_true)) ' and sampled g: \mu=' num2str(mean(g.d)) ' | \sigma='   num2str(std(g.d))]);

plot(g.x,g.y,'or')
for i_sim=2:m*k
    subplot(k,m,i_sim); pcolor(grid{end}.x,grid{end}.y,gG{end}.m{i_sim-1}); xlabel('x[m]');ylabel('y[m]');
    title(['Simulated Electrical Conductivity gG: \mu=' num2str(mean2(gG{end}.m{i_sim-1})) ' | \sigma=' num2str(std2(gG{end}.m{i_sim-1})) ]);
    shading flat;colorbar;caxis([min(g_true(:)) max(g_true(:))]);axis tight
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1); hist(g_true(:));title('True Electrical Conductivity g_{true}')
subplot(1,3,2); hist(g.d(:));title('Sampled True Electrical Conductivity g')
subplot(1,3,3); hist(gG{end}.m{i_sim-1}(:));title('Reproduced Electrical Conductivity gG')



figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1); pcolor(grid{end}.x,grid{end}.y,K_true);xlabel('x[m]');ylabel('y[m]');title('True Hydraulic conductivity K_{true}');shading flat;colorbar;
subplot(2,1,2); pcolor(grid{end}.x,grid{end}.y,kK{end}.m);xlabel('x[m]');ylabel('y[m]');title('Simulated Hydraulic conductivity kK');shading flat;colorbar;


%% Saveit