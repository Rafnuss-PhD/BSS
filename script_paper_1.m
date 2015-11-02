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
parm.scale      = 10;
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
xlabel('Grid size'); ylabel('Time [s]')
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
parm.n_realisation  = 2; parm.name = 'kcst1-2';
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
parm.n_realisation  = 30; parm.name = 'kcst0-30';
BSGS(sigma,Sigma,sigma_true,grid,parm);

figure; hold on;
list = ls('result/SimilarToPaolo/kcst*');
for i=1:size(list,1)
    load(strtrim(list(i,:)))
    gg=[grid{parm.scale}];
    if parm.cstk
        LineSpec = '--x';
    else
        LineSpec = '-o';
    end
    plot([gg.nxy],cell2mat(t.scale),LineSpec,'DisplayName',['Same path: ' num2str(parm.cstk) '. Number of realisation: ', num2str(parm.n_realisation)])
end
xlabel('Grid size'); ylabel('Time [s]')
ax = gca; ax.XTick=sort([ax.XTick,[gg.nxy]]);
legend(ax,'show')
ax.Yscale='log';
ax.XScale='log';axis tight;

figure; hold on; xlabel('Number of Realisation'); ylabel('Time [min]')
%ax = gca; ax.YScale='log'; 
list = ls('result/SimilarToPaolo/kcst*');
data_0=[];data_1=[];
for i=1:size(list,1)
    load(strtrim(list(i,:)))
    if parm.cstk
        LineSpec = 'x';
        data_1=[data_1;parm.n_realisation,t.global];
    else
        LineSpec = 'o';
        data_0=[data_0;parm.n_realisation,t.global];
    end
    scatter(parm.n_realisation,t.global/60,LineSpec,'LineWidth',2)
end

fit_0=fit(data_0(:,1),data_0(:,2),'a*x+b','StartPoint',[50 100]);
plot(fit_0)
fit_1=fit(data_1(:,1),data_1(:,2),'a*x+b','StartPoint',[50 100]);
plot(fit_1)



%% Multi-grid scale to optimize time vs exploring space.
% comnpute from previous point the number of realisation which last 10min
load('result/SimilarToPaolo/SimilarToPaolo_Zstd_corrected.mat');parm.gen = gen;
parm.neigh = 1; 
parm.nb_neigh  = [0 0 0 0 0; 5 5 5 5 5];
parm.scale = 1:8;
parm.cstk = true;
parm.n_realisation  = 20; parm.name = 'kcst1-11';
BSGS(sigma,Sigma,sigma_true,grid,parm);
parm.cstk = false;
parm.n_realisation  = 4; parm.name = 'kcst0-4';
BSGS(sigma,Sigma,sigma_true,grid,parm);

fieldname = {'kcst0-4_2015-11-02_12-15-06', 'kcst1-11_2015-11-02_11-54-02','kcst8-9_2015-11-02_14-29-51'};
figure; hold on; axis equal
for i= 1:numel(fieldname)
    load(['result/Cstk-space-exploration/' fieldname{i}])
    data=[];
    for i_realisation =1:parm.n_realisation
        data= [data, Y{end}.m{i_realisation}(:)];
    end
    D = pdist(data');
    D_y = cmdscale(D);
    scatter(D_y(:,1),D_y(:,2),'DisplayName',['Identical path : ', num2str(parm.cstk), ' in ' num2str(round(t.global/60)) 'min'])
end
ax = gca;legend(ax,'show')


figure; 
sum_field=cell(3,2);
for i= 1:numel(fieldname)
    load(['result/Cstk-space-exploration/' fieldname{i}])
    sum_field{i,1}=nan(grid{end}.ny,grid{end}.nx);
    sum_field{i,2}=nan(grid{end}.ny,grid{end}.nx);
    for x = 1:grid{end}.nx
        for y = 1:grid{end}.ny
            data=[];
            for i_realisation =1:parm.n_realisation
                data= [data, Y{end}.m{i_realisation}(y,x)];
            end
            sum_field{i,1}(y,x) = mean(data);
            sum_field{i,2}(y,x) = std(data);
        end
    end
    
end

figure;
c_axis_1=[min(min([sum_field{:,1}])) max(max([sum_field{:,1}]))];
c_axis_2=[min(min([sum_field{:,2}])) max(max([sum_field{:,2}]))];
subplot(3,2,1);imagesc(sum_field{1,1}); caxis(c_axis_1)
subplot(3,2,2);imagesc(sum_field{1,2}); caxis(c_axis_2)
subplot(3,2,3);imagesc(sum_field{2,1}); caxis(c_axis_1)
subplot(3,2,4);imagesc(sum_field{2,2}); caxis(c_axis_2)
subplot(3,2,5);imagesc(sum_field{3,1}); caxis(c_axis_1)
subplot(3,2,6);imagesc(sum_field{3,2}); caxis(c_axis_2)


%% Multi-grid scale to optimize time vs exploring space.
parm.scale = 1:6;
n_realisation = 1:10;
cstk_s        = 0:7;
[N_realisation, cstk_S]=meshgrid(n_realisation,cstk_s);
parfor i=1:numel(N_realisation)
    parm1=parm;
    parm1.cstk_s=cstk_S(i);
    parm1.n_realisation  = N_realisation(i);%n_realisation(i);
    parm1.name = ['kcst' num2str(parm1.cstk_s) '-' num2str(parm1.n_realisation)];
    BSGS(sigma,Sigma,sigma_true,grid,parm1);
end



%ax = gca; ax.YScale='log'; 
list = ls('result/Cstk-space-exploration-2/kcst*');
data=[];
for i=1:size(list,1)
    load(['result/Cstk-space-exploration-2/' strtrim(list(i,:))])
    data=[data;parm.n_realisation,parm.cstk_s,t.global/60, std2([Y{end}.m{:}])];
end

figure; hold on; xlabel('Number of Realisation'); ylabel('Scale from which path is constant'); title('Time'); hold on; colorbar
F = scatteredInterpolant(data(:,1),data(:,2),data(:,3));
[xq,yq] = meshgrid(min(data(:,1)):max(data(:,1)), min(data(:,2)):max(data(:,2)));
contourf(xq,yq,F(xq,yq));
scatter(data(:,1),data(:,2),[],data(:,3),'MarkerEdge','k')



%% Neighbouring improvement
addpath(genpath('./.')); load('result/Random_Neigh/Random_10x10_2015-11-02_17-00.mat')
parm.gen            = gen;
parm.nscore         = true;
parm.likelihood     = false;
parm.n_realisation  = 1;

parm.scale          = 1:6;
parm.neigh          = true;
parm.nb_neigh       = [0 0 0 0 0; 5 5 5 5 5];
%parm.k.sb.nx        = ceil(grid{end}.x(end)/gen.covar.modele(1,2)*3);
%parm.k.sb.ny        = ceil(grid{end}.y(end)/gen.covar.modele(1,3)*3);
parm.scale          = 10;
parm.gen.covar.modele  = [4 40 4 0;1 1 1 0];
parm.gen.covar.c  = [.99;.01];

