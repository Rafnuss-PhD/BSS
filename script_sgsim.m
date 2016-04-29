clear all;clc;

% General setting
grid_gen.sx = 7; % (65 x 65)
grid_gen.sy = 7;
grid_gen.x = 0:2^grid_gen.sx;
grid_gen.y = 0:2^grid_gen.sy;
sigma_true=[]; % no hard data yet
sigma.d=[];
sigma.x=[];
sigma.y=[];
sigma.n=0;
parm.gen.covar.modele = [4 20 20 0]; % gaussian anisotrope 20x10, no nuggets
parm.gen.covar.c = 1;
parm.saveit = false;
parm.nscore = 0;
parm.par = 0;
parm.n_realisation  = 1;
parm.cstk = true;
parm.seed = 1000;
parm.varcovar = 1;
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.neigh = 1;

% the true covariance matrix for unconditional model
X = repmat(grid_gen.x,numel(grid_gen.y),1);
Y = repmat(grid_gen.y',1,numel(grid_gen.x));
profile resume;
[Res, ~, ~, parm, ~] = SGSIM(sigma,sigma_true,grid_gen,parm); % fake one to generate cx
profile off; profile viewer
CY=covardm_perso([X(:) Y(:)],[X(:) Y(:)],parm.gen.covar.modele,parm.gen.covar.c,parm.k.cx);
nu = @(res,res_ref) full( sqrt(sum((res(:)-res_ref(:)).^2)) / sqrt( sum(res_ref(:).^2)));


%% Test 1 : How many neighbors point do we need for the difference path method
parm.k.wradius = Inf; 
parm.saveit = true;
nb_neigh=[3 5 10 30];
path = {'row-by-row','spiralout','spiralin','maximize','random'};

parfor i_nb_neigh=1:numel(nb_neigh)
    for i_path=1:numel(path)
        parm_t = parm;
        parm_t.path=path(i_path);
        parm_t.k.nb_neigh = [0 0 0 0 0; repmat(nb_neigh(i_nb_neigh),1,5)];
        parm_t.name = path+ '-'+ num2str(parm_t.k.nb_neigh);
        SGSIM(sigma,sigma_true,grid_gen,parm_t);
    end
end

for i_path=1:numel(Res)
    nu_nb_neigh(i_path) = nu(Res{i_path}{end}.varcovar.CY,CY);
    t_nb_neigh(i_path) = t{i_path}.global;
    var_v_nb_neigh(i_path,:) = Res{i_path}{end}.varcovar.vario_1d_v;
    var_h_nb_neigh(i_path,:) = Res{i_path}{end}.varcovar.vario_1d_h;
end


%% Test 2 neighborhod

nb_neigh=[3 5 10 30 100];
for i_nb_neigh=1:numel(nb_neigh)
    parm.k.nb_neigh = [0 0 0 0 0; repmat(nb_neigh(i_nb_neigh),1,5)];
    [Res{i_nb_neigh}, t{i_nb_neigh}, ~, parm, ~] = SGSIM(sigma,sigma_true,grid_gen,parm);
end

for i_nb_neigh=1:numel(Res)
    nu_nb_neigh(i_nb_neigh) = nu(Res{i_nb_neigh}{end}.varcovar.CY,CY);
    t_nb_neigh(i_nb_neigh) = t{i_nb_neigh}.global;
    var_v_nb_neigh(i_nb_neigh,:) = Res{i_nb_neigh}{end}.varcovar.vario_1d_v;
    var_h_nb_neigh(i_nb_neigh,:) = Res{i_nb_neigh}{end}.varcovar.vario_1d_h;
end

figure(1); hold on;
plot(nb_neigh,nu_nb_neigh*100,'-x');ylabel('\nu [%]'); yyaxis right
plot(nb_neigh,t_nb_neigh,'-x'); ylabel('Time [sec]')
xlabel('Number of neighbors/quadrant')

figure(2); hold on;
myvario = @(h,range) semivariogram1D(h,parm.k.var,range,parm.k.model(1),0);
plot(grid_gen.x, myvario(grid_gen.x, parm.k.model(3)),'--k')
for i_nb_neigh=1:numel(Res)
    plot(grid_gen.x, parm.k.var-Res{i_nb_neigh}{end}.varcovar.vario_1d_h.^2)
end
plot(grid_gen.x, parm.k.var-Res_ref{end}.varcovar.vario_1d_h)
legend({'theorical', mat2cell(nb_neigh,1), 'Inf' })

figure(3); hold on;
subplot(3,2,1); imagesc(Res_ref{end}.varcovar.CY);
for i_nb_neigh=1:numel(Res)
   subplot(3,2,i_nb_neigh+1); imagesc(Res{i_nb_neigh}{end}.varcovar.CY-Res_ref{end}.varcovar.CY); caxis([-.1, .1])
   legend(nb_neigh(i_nb_neigh)); colorbar;
end






%% Test 2: path 

parm.gen.covar.modele = [4 round(grid_gen.x(end)/3) round(grid_gen.x(end)/3)/2 0];
parm.gen.covar.modele = [4 5 1 0];

parm.k.nb_neigh = [0 0 0 0 0; Inf Inf Inf Inf Inf]; 


[Res_ref, ~, ~, parm, ~] = SGSIM(sigma,sigma_true,grid_gen,parm);
figure(10); imagesc(Res_ref{end}.varcovar.CY);

figure(11); imagesc(Res_ref{end}.varcovar.vario_2d);

