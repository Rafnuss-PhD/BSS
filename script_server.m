clear all; close all
load('result-BSGS/GEN-Run_1_2017-05-07_14-37');
[Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
parm.k.covar = gen.covar;
parm.kernel=kern;

parm.unit='';
parm.nscore = 1;
parm.cstk = true;
parm.seed = 'shuffle';
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'minkmex'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 0;
parm.plot.krig=0;

parm.k.n = 40;

parm.n_realisation  = 0;

parm.par = 0;
parm.par_n = 4;

% use the log of hyd. cond.
Klog=K;
Klog.d=log(Klog.d);

%% work on the weight
% parm.aggr.method='cst';
% parm.aggr.T = .1:.1:.9;

parm.aggr.method='step';
parm.aggr.T = .006;

% parm.aggr.method='linear';
% parm.aggr.T = [ .1 .9];
%             
% parm.aggr.method='sigmoid';
% parm.aggr.T = [ .1 .9];


parm.aggr.sum = 1;
parm.n_realisation  = parm.par_n*numel(parm.aggr.T)*8;


disp('Setup ok, Start')
[Res,~,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);
disp('Run finished')

% Sigma_d = Sigma.d(any( bsxfun(@eq,Sigma.X(:),Res.X(:)') & bsxfun(@eq,Sigma.Y(:),Res.Y(:)' ),2));
id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));

E1 = nan(sum(id),parm.n_realisation);
E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
    
XY = kern.XY;
Sigma_d = Sigma.d(:);
Res_m_ns = Res.m_ns;
Res_m = Res.m;
dens = parm.kernel.dens(:);
parfor i_realisation=1:parm.n_realisation
   gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
   E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
   E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
end

try
    parfor i_x = 1:size(parm.aggr.T,1)
        idx=i_x:numel(parm.aggr.x):parm.n_realisation;
        OF1(i_x) = sqrt(mean(mean(E1(:,idx),2).^2));
        OF2(i_x) = sqrt(mean(mean(E2(:,idx),2).^2));
    end
catch
    save(['result-BSS/' filename],'Res','parm','kern','Nscore','E1','E2')
end

save(['result-BSS/' filename],'Res','parm','kern','Nscore','E1','E2','OF1','OF2')
disp('file written')
%save('result-BSS/ResAB_A_1.mat','parm','E1','E2')
% save('C:\Users\Raphael\Desktop\to be compress\ResRAD_AB.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')

