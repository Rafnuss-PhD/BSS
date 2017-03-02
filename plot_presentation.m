test = [255, 255, 204;161, 218, 180;65, 182, 196;44, 127, 184;8, 104, 172;37, 52, 148];
test = [interp1(linspace(1,64,size(test,1)), test(:,1), 1:64)' interp1(linspace(1,64,size(test,1)), test(:,2), 1:64)' interp1(linspace(1,64,size(test,1)), test(:,3), 1:64)'];


%% FRONT PAGE
grid_gen.sx = 8; % (65 x 65)
grid_gen.sy = 8;
grid_gen.x = 0:2^grid_gen.sx;
grid_gen.y = 0:2^grid_gen.sy;

% Default variogram
parm.gen.covar(1).model = 'spherical';
parm.gen.covar(1).range = [100 100];
parm.gen.covar(1).azimuth = 0;
parm.gen.covar(1).c0 = .98;
parm.gen.covar(1).alpha=1;
% parm.gen.covar(2).model = 'nugget';
% parm.gen.covar(2).range = [15 15];
% parm.gen.covar(2).azimuth = 0;
% parm.gen.covar(2).c0 = 0.02;

addpath(genpath('./.'))
parm = kriginginitiaite(parm);
field = fftma_perso(parm.k.covar, grid_gen);

imagesc(field)
colormap(test/255)


%% SGSIM
grid_gen.sx = 4; % (65 x 65)
grid_gen.sy = 4;
grid_gen.x = 0:2^grid_gen.sx;
grid_gen.y = 0:2^grid_gen.sy;
[grid_gen.X, grid_gen.Y] = meshgrid(grid_gen.x, grid_gen.y);
X.d=[];X.n=0;X.x=[];X.y=[];
parm.saveit = false;
parm.nscore = 0;
parm.par = 0;
parm.n_realisation  = 1;
parm.cstk = true;
parm.seed = 1111;
parm.varcovar = 1;
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.neigh = 0;
parm.k.wradius=Inf;
parm.k.nb_neigh=[0 0 0 0 0; 4 0 0 0 0];

% Default variogram
parm.gen.covar(1).model = 'spherical';
parm.gen.covar(1).range = [5 5];
parm.gen.covar(1).azimuth = 0;
parm.gen.covar(1).c0 = .98;
parm.gen.covar(1).alpha=1;
% parm.gen.covar(2).model = 'nugget';
% parm.gen.covar(2).range = [15 15];
% parm.gen.covar(2).azimuth = 0;
% parm.gen.covar(2).c0 = 0.02;

parm = kriginginitiaite(parm);

parm.path_random=1;
parm.path = 'maximize';
parm.plot.krig=0;
parm.plot.presentation=1;

[ResInf, t, k, parm, filename] = SGSIM(X,grid_gen,parm);

%%

[h, w,p] = size(ResInf{end}.F(1).cdata);  % use 1st frame to get dimensions
hf = figure; 
% resize figure based on frame's w x h, and place at (150, 150)
set(hf, 'position', [150 150 w h]);
axis off
% movie(hf,Res{end}.F,3,2);
%mplay(f)
myVideo = VideoWriter('SGSIM_MG_4nodes.avi','Uncompressed AVI');
myVideo.FrameRate = 3;
open(myVideo); writeVideo(myVideo,ResInf{end}.F); close(myVideo);

%%

figure(2);clf
[xInf,yInf]=variogram_gridded_perso(ResInf{end}.m{end});
[x8,y8]=variogram_gridded_perso(Res8{end}.m{end});
subplot(2,1,1); hold on;
plot(grid_gen.x, 1-parm.k.covar.g(grid_gen.x/parm.k.covar.range(1)),'--k','linewidth',2);
plot(grid_gen.x, xInf,'linewidth',2); 
plot(grid_gen.x, x8,'linewidth',2); 
ylabel('Vertical variogram'); xlabel('Lag-distance')
set(gca,'color','none')
set(gca,'Linewidth',1)
xlim([0 10])

subplot(2,1,2);hold on;
plot(grid_gen.y, 1-parm.k.covar.g(grid_gen.y/parm.k.covar.range(2)),'--k','linewidth',2);
plot(grid_gen.y,yInf,'linewidth',2);
plot(grid_gen.y,y8,'linewidth',2);
ylabel('Horizontal variogram'); xlabel('Lag-distance'); 
set(gca,'color','none')
set(gca,'Linewidth',1)
xlim([0 10])


legend('Theorical Variogram','Limited Neighborhood','Global Neighborhood')

export_fig test.png -transparent


%%
figure(1);clf; hold on;
plot(0:5,parm.k.covar.g(0:0.2:1.1),'o-k')
set(gca,'Ytick',sort(parm.k.covar.g(0:0.2:1.1)));
set(gca,'Yticklabel',{'C_5','C_4','C_3','C_2','C_1','C_0'});
xlabel('Distance'); ylabel('Covariogram')



