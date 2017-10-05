clear all; close all; addpath(genpath('./.'));
load('result-BSGS/GEN-Run_1_2017-05-07_14-37');
[Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
parm.k.covar = gen.covar;
parm.kernel=kern;

parm.unit='';
parm.nscore = 1;
parm.cstk = true;
parm.seed = 'shuffle';
% parm.scale=[1:grid_gen.sx;1:grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 0;
parm.plot.krig=0;


% parm.k.n = 40;
parm.par = 1;
parm.par_n = 4;

% use the log of hyd. cond.
Klog=K;
Klog.d=log(Klog.d);

%% work on the weight
% parm.aggr.method='cst';
% parm.aggr.T = .5;

filename='ResStepVarious';
parm.aggr.method='step';
parm.aggr.T = [0  .001 .002 .003 .06 .1 .5 1]';

% parm.aggr.method='linear';
% parm.aggr.T = [ .1 .9];
    
% filename='ResSigmoid06';
% parm.aggr.method='sigmoid';
% parm.aggr.T = [ .06 Inf ; .06 1000; .06  100; .06  50; .06  20; .06  10];


parm.aggr.sum = 1;
parm.n_realisation  = parm.par_n*size(parm.aggr.T,1);


disp('Setup ok, Start')
[Res,~,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);
disp('Run finished')

% Sigma_d = Sigma.d(any( bsxfun(@eq,Sigma.X(:),Res.X(:)') & bsxfun(@eq,Sigma.Y(:),Res.Y(:)' ),2));
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));

E1 = nan(sum(id),parm.n_realisation);
E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
    
XY = kern.XY;
Sigma_d = Sigma.d(:);
Res_m_ns = Res.m_ns;
Res_m = Res.m;
dens = parm.kernel.dens(:);
for i_realisation=1:parm.n_realisation
   gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
   E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
   E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
end
disp('Error Computed')


try
    clear OF1 OF2
    for i_ab = 1:size(parm.aggr.T,1)
        i_t=i_ab:size(parm.aggr.T,1):parm.n_realisation;
        OF1(i_ab) = sqrt(mean(mean(E1(:,i_t),2).^2));
        OF2(i_ab) = sqrt(mean(mean(E2(:,i_t),2).^2));
    end
catch
    %save(['result-BSGS/' filename],'Res','parm','kern','Nscore','E1','E2')
end

%save(['result-BSGS/' filename],'Res','parm','kern','Nscore','E1','E2','OF1','OF2')
disp('file written')

figure(1); clf; hold on; 
scatter(OF1,OF2,[],1:numel(OF1),'filled')

filename = {'ResCst11','ResA0','ResA1','ResCst5','ResCst111','ResCst5'};%
for i=1:numel(filename)-1
    load(['result-BSGS\' filename{i} '.mat'], 'OF1', 'OF2')
    scatter(OF1,OF2,'filled','r')
end
