clear all;clc;

% General setting
grid_gen.sx = 6; % (65 x 65)
grid_gen.sy = 6;
grid_gen.x = 0:2^grid_gen.sx;
grid_gen.y = 0:2^grid_gen.sy;
sigma_true=[]; % no hard data yet
sigma.d=[];
sigma.x=[];
sigma.y=[];
sigma.n=0;
parm.gen.covar(1).model = 'spherical';
parm.gen.covar(1).range = [15 15];
parm.gen.covar(1).azimuth = 0;
parm.gen.covar(1).c0 = 1;
parm.saveit = false;
parm.nscore = 0;
parm.par = 0;
parm.n_realisation  = 1;
parm.cstk = true;
parm.varcovar = 1;
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.neigh = 0;
parm.seed = 'shuffle';

parm = kriginginitiaite(parm);


%%

col={[0 113 188]/255, [216 82 24]/255, [236 176 31]/255 [125 46 141]/255};

% path
path = {'linear''linear','maximize'};
path_random = [0 1 1];
parm.k.wradius=Inf;
vario = {'exponential', 'gaussian','spherical', 'stable', 'hyperbolic','k-bessel', 'cardinal sine'};
neigh = [20];

% Default variogram
parm.gen.covar(1).model = '';
parm.gen.covar(1).range = [15 15];
parm.gen.covar(1).azimuth = 0;
parm.gen.covar(1).c0 = .98;
parm.gen.covar(2).model = 'nugget';
parm.gen.covar(2).range = [ 0 0];
parm.gen.covar(2).azimuth = 0;
parm.gen.covar(2).c0 = 0.02;
parm.gen.covar(1).alpha = 1;

% the true covariance matrix for unconditional model
[X1, Y1] = meshgrid(grid_gen.x,grid_gen.y);
X = repmat(0:(grid_gen.x(end)*2),numel(grid_gen.y)*2-1,1);
Y = repmat((0:grid_gen.y(end)*2)',1,numel(grid_gen.x)*2-1);

for j=1:numel(vario)
    parm.gen.covar(1).model = vario{j};
    parm = kriginginitiaite(parm);
    CY_true{j} = covardm_perso([X1(:) Y1(:)],[X1(:) Y1(:)],parm.k.covar);
    vario_true{j} = reshape(covardm_perso([X(:) Y(:)],[grid_gen.x(end) grid_gen.y(end)],parm.k.covar),numel(grid_gen.y)*2-1,numel(grid_gen.x)*2-1);
    CY_true_ss{j} = sqrt( sum(CY_true{j}(:).^2));
end
clear X Y X1 Y1

err_frob_fx = @(CY,j) full( sqrt(sum((CY(:)-CY_true{j}(:)).^2)) / CY_true_ss{j});


% Simulation
for w = 1:numel(path)
    for j=1:numel(vario)
        parm.path = path{w};
        parm.gen.covar(1).model=vario{j};
        parm.path_random = path_random(w);
        parm = kriginginitiaite(parm);
        for i= 1:numel(neigh)
            parm.k.nb_neigh = [0 0 0 0 0; neigh(i) 0 0 0 0];
            [Res,~, ~, ~, ~] = SGSIM(sigma,sigma_true,grid_gen,parm);
            [CY{w,i,j},~,~,~,vario_sim{w,i,j}] = varcovar(Res, 'vario', 1);
        end
    end
end


% Error
for w = 1:numel(path)
    for i= 1:numel(neigh)
       err{j}(w,i) = err_frob_fx(CY{w}{i});
    end
end



% Figure Error
figure(2);clf;
 for j=1:numel(vario_type)
     subplot(2,2,j); hold on
     for w = 1:numel(path)
         plot(neigh, reshape(err(j,w,:),numel(neigh),1))
     end
     xlabel(sprintf('Neighborhood size \n %s variogram', vario_name{j}));
     ylabel('Frobenium Norm Error')
     legend(path)
 end

% Figure Vario
figure(3);clf;
 neigh_s = strread(num2str(neigh),'%s');
 neigh_s = {'true' neigh_s{:}};
 path_s={'Row-by-Row' 'Random' 'Multigrid'};
for w=1:numel(path)
    for j=1:numel(vario)
        subplot(numel(path),numel(vario),j+ (w-1)*numel(vario)); hold on;
        plot(vario_true{j}(end/2+.5,end/2+.5:end),'--k')
%         for i=1:numel(neigh)
%             plot(vario_sim{w,i,j}{1}.h)
%         end
        %legend(neigh_s)
        axis tight
        if (w==numel(path))
            xlabel(['Distence\newline' vario{j}]);
        end
        if (j==1)
            ylabel(['Covariogramm\newline' path_s(w)]);
        end
    end
end
