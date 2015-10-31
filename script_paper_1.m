% Script paper 1

%% Plot for Smart neighbouring: 
% This simulation are meant to shows that for large grid, smart
% neighbooring largly improve simulations computational cost


% load the input data
addpath(genpath('./.')) 
load('result/Neigh/Very_large_1millionscells_2015-10-27_10-52.mat')

% Define std parameter
parm.n_realisation  = 1;
parm.scale          = 6;

parm.neigh          = true;
parm.cstk           = true;
parm.nscore         = true;

% Saving
parm.saveit         = false;
parm.gen            = gen;

% Parameter : 
% - nb of kriging point: min,max per quadrant | sum of max
% - size of superblock grid
% - size of windows search.
parm.nb_neigh   = [0 0 0 0; 0 0 0 0];
parm.k.sb.nx    = ceil(grid{end}.x(end)/gen.covar.modele(1,2)*3);
parm.k.sb.ny    = ceil(grid{end}.y(end)/gen.covar.modele(1,3)*3);
parm.scale      = 6;
parm.gen.covar.modele  = [4 40 4 0;1 1 1 0];
parm.gen.covar.c  = [.99;.01];


% Varying nb of kriging point
nb_neigh=[4 8 15];
nb_neigh_hd=[10 10 10];
neigh = [0 1];
for i=1:numel(neigh)
    parm.neigh = neigh(i);
    name_b = ['neigh' num2str(neigh(i))];
    for u = 1:numel(nb_neigh)
        parm.name = [name_b, '_nbk' sprintf('%02d',nb_neigh(u))];
        parm.nb_neigh(2,:)   = [repmat(nb_neigh(u),1,4) nb_neigh_hd(u)];
        BSGS(sigma,Sigma,sigma_true,grid,parm);
    end
end

figure; hold on;
list = ls('result/Neigh/neigh*');
for i=1:size(list,1)
    load(strtrim(list(i,:)))
    gg=[grid{parm.scale}];
    if parm.neigh
        LineSpec = '--x';
    else
        LineSpec = '-o';
    end
    plot([gg.nxy],cell2mat(t.scale),LineSpec,'DisplayName',['Number of kriging points: ', num2str(sum(parm.nb_neigh(2,:)))])
end
xlabel('Grid size'); label('Time [s]')
ax = gca; ax.XTick=sort([ax.XTick,[gg.nxy]]);
legend(ax,'show')



%% Similar to Paolo
load('result/SimilarToPaolo/SimilarToPaolo_Zstd_corrected.mat');parm.gen = gen;

% With-without Neighbouring / Multi-grid
parm.nb_neigh  = [0 0 0 0 0; 5 5 5 5 5];
parm.neigh = 0; parm.scale = 10; parm.name = 'MG0-neigh0';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.neigh = 0; parm.scale = 1:10; parm.name = 'MG1-neigh0';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.neigh = 1; parm.scale = 10; parm.name = 'MG0-neigh1';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.neigh = 1; parm.scale = 1:10; parm.name = 'MG1-neigh1';
BSGS(sigma,Sigma,sigma_true,grid,parm);

result=[382 211 318 149]; % in secondes.

% Constent Path
parm.cstk = true;
parm.neigh = 1; 
parm.scale = 1:10;
parm.n_realisation  = 1; parm.name = 'kcst1-1';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.n_realisation  = 5; parm.name = 'kcst1-5';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.n_realisation  = 20; parm.name = 'kcst1-20';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.n_realisation  = 100; parm.name = 'kcst1-100';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.cstk = false;
parm.n_realisation  = 1; parm.name = 'kcst0-1';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.n_realisation  = 2; parm.name = 'kcst0-2';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.n_realisation  = 10; parm.name = 'kcst0-10';
BSGS(sigma,Sigma,sigma_true,grid,parm);