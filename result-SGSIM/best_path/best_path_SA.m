clear all;clc;
addpath(genpath('./.'))

% General setting
grid_gen.sx = 6; % (65 x 65)
grid_gen.sy = 6;
grid_gen.x = 0:2^grid_gen.sx;
grid_gen.y = 0:2^grid_gen.sy;
sigma.d=[];sigma.x=[];sigma.y=[];sigma.n=0;
parm.saveit = false;
parm.nscore = 0;
parm.par = 0;
parm.n_realisation  = 1;
parm.cstk = true;
parm.varcovar = 1;
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.neigh = 0;
parm.seed = 'shuffle';
parm.k.wradius=Inf;

% Default variogram
parm.gen.covar(1).model = 'spherical';
parm.gen.covar(1).range = [15 15];
parm.gen.covar(1).azimuth = 0;
parm.gen.covar(1).c0 = .98;
parm.gen.covar(2).model = 'nugget';
parm.gen.covar(2).range = [15 15];
parm.gen.covar(2).azimuth = 0;
parm.gen.covar(2).c0 = 0.02;
parm.gen.covar(1).alpha = 1;

parm.k.nb_neigh = [0 0 0 0 0; 10 0 0 0 0];
parm = kriginginitiaite(parm);

% the true covariance matrix for unconditional model
[grid_gen.X, grid_gen.Y] = meshgrid(grid_gen.x,grid_gen.y);
X = repmat(0:(grid_gen.x(end)*2),numel(grid_gen.y)*2-1,1);
Y = repmat((0:grid_gen.y(end)*2)',1,numel(grid_gen.x)*2-1);

CY_true = covardm_perso([grid_gen.X(:) grid_gen.Y(:)],[grid_gen.X(:) grid_gen.Y(:)],parm.k.covar);
vario_true = reshape(covardm_perso([X(:) Y(:)],[grid_gen.x(end) grid_gen.y(end)],parm.k.covar),numel(grid_gen.y)*2-1,numel(grid_gen.x)*2-1);
CY_true_ss = sqrt( sum(CY_true(:).^2));
err_frob_fx = @(CY) full( sqrt(sum((CY(:)-CY_true(:)).^2)) / CY_true_ss);
clear X Y

p_rbr = 1:numel(grid_gen.X);
parm.path='maximize';
parm.path_random=1;
[Res,~, ~, ~, ~] = SGSIM(sigma,grid_gen,parm);
p_mg = Res{1}.path;

%% Version 

fun = @(path) err_frob_fx(varcovar(SGSIM(sigma,grid_gen,parm_path(parm,path)), 'vario', 0));

options = optimoptions(@simulannealbnd, ...
    'AnnealingFcn',@AnnealingFcn,...
    'MaxIterations',10000,...
    'FunctionTolerance',0,...
    'MaxTime',60*2,...
    'MaxStallIterations', 500),...
    'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
    
parfor i=1:4
    p_r(i,:)=randperm(numel(grid_gen.X));
    f_r(i) = fun(p_r(i,:))
    [p_sr(i,:),f_sr(i)] = simulannealbnd(fun,p_r(i,:),0,2^32,options);
    %[p_smg(i,:),f_smg(i)] = simulannealbnd(fun,p_mg,0,2^32,options);
end

f_rbr = fun(p_rbr);
f_mg = fun(p_mg);

%save('best_path_SA.m','p_r','p_sr','p_smg','p_mg','p_rbr','f_r','f_sr','f_smg','f_rbr','f_mg')

%% Figure 

% Histogram
figure(3); clf; hold on
%histogram(fval_p0)
histogram(fval)
plot(fval_rbr,0,'o')
plot(fval_mg,0,'o')
legend('Random path','Random converged','Row-by-row','Multi-grid')
xlabel('Frobenium error')


% get distribution
for i=1:20
    Res=SGSIM(sigma,grid_gen,parm_path(parm,p(i,:)));
    dist{i} = Res{1}.dist;
end
Res_rbr=SGSIM(sigma,grid_gen,parm_path(parm,p_rbr));
Res_mg=SGSIM(sigma,grid_gen,parm_path(parm,p_mg));
figure; hold on;
for j=1:20
    for i = 1:81
        s(i) = mean(dist{j}(dist{j}(:,2)==i,1));
    end
    plot(s);
end
for i = 1:81
    s1(i) = mean(Res_rbr{1}.dist(Res_rbr{1}.dist(:,2)==i,1));
    s2(i) = mean(Res_mg{1}.dist(Res_mg{1}.dist(:,2)==i,1));
end
plot(s1);plot(s2);


% Grid
m=nan(Res{end}.nx,Res{end}.ny);  i=1;
x=p_mg;
for i_scale=1:1
    for i_pt = 1:Res{i_scale}.sim.n
        m(Res{i_scale}.varcovar_id(x(i_pt))) = i_pt;
    end
end
figure(35); clf; hold on;
imagesc(grid_gen.x,grid_gen.y,m,'AlphaData',~isnan(m)); hold on; plot(sigma.x+1,sigma.y+1,'sr','MarkerFaceColor','r','MarkerSize',20); 
axis tight equal; xlabel('x');ylabel('y'); axis tight; 

% nearest node
for i_p = 1:12
    [s_mean_r(i_p) s_min_r(i_p), s_un_r(i_p)] = nearestnode(p_r(i_p,:),parm);
    [s_mean_sr(i_p) s_min_sr(i_p), s_un_sr(i_p)] = nearestnode(p_sr(i_p,:),parm);
    [s_mean_smg(i_p) s_min_smg(i_p), s_un_smg(i_p)] = nearestnode(p_smg(i_p,:),parm);
end
[s_mean_rbr s_min_rbr, s_un_rbr] = nearestnode(p_rbr,parm);
[s_mean_mg s_min_mg, s_un_mg] = nearestnode(p_mg,parm);

figure(45);clf;
subplot(2,2,1); hold on;
scatter(s_min_r,f_r)
scatter(s_min_sr,f_sr)
scatter(s_min_smg,f_smg)
plot(s_min_rbr,f_rbr,'o')
plot(s_min_mg,f_mg,'o')
xlabel('Average nearest node');ylabel('Frebenius Norm')
legend('Random','Random converged','MG converged','Row-by-row','Multi-grid')
subplot(2,2,2); hold on;
scatter(s_mean_r,f_r)
scatter(s_mean_sr,f_sr)
scatter(s_mean_smg,f_smg)
plot(s_mean_rbr,f_rbr,'o')
plot(s_mean_mg,f_mg,'o')
legend('Random','Random converged','MG converged','Row-by-row','Multi-grid')
xlabel('Average Mean of distence to neighboors');ylabel('Frebenius Norm')
subplot(2,2,3); hold on;
scatter(s_un_r,f_r)
scatter(s_un_sr,f_sr)
scatter(s_un_smg,f_smg)
plot(s_un_rbr,f_rbr,'o')
plot(s_un_mg,f_mg,'o')
legend('Random','Random converged','MG converged','Row-by-row','Multi-grid')


% kriging variance
for i_p = 1:12
    Res=SGSIM(sigma,grid_gen,parm_path(parm,p_r(i_p,:)));
    s_s_r(i_p,:) = Res{end}.s;
    Res=SGSIM(sigma,grid_gen,parm_path(parm,p_sr(i_p,:)));
    s_s_sr(i_p,:) = Res{end}.s;
    Res=SGSIM(sigma,grid_gen,parm_path(parm,p_smg(i_p,:)));
    s_s_smg(i_p,:) = Res{end}.s;
end
Res=SGSIM(sigma,grid_gen,parm_path(parm,p_rbr));
s_s_rbr = Res{end}.s;
Res=SGSIM(sigma,grid_gen,parm_path(parm,p_mg));
s_s_mg = Res{end}.s;

figure(45);clf; hold on;
plot(mean(s_s_r))
plot(mean(s_s_sr))
plot(mean(s_s_smg))
plot(s_s_rbr)
plot(s_s_mg)