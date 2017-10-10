<<<<<<< HEAD
addpath(genpath('./.'));
load('result-BSS/GEN-Run_1_2017-05-07_14-37');

Nscore = nscore(kern, struct('nscore', 1), 0); %Prim.d, kern.axis_prim, 'pchip', 'pchip', parm.plot.ns

[~,s_id]=min(bsxfun(@minus,kern.axis_sec,Sigma.d(:)).^2,[],2);
sec_pdf = kern.dens(:,s_id);
sec.pdf = bsxfun(@times, sec_pdf, 1./sum(sec_pdf));
sec.axis = Nscore.forward(kern.axis_prim);

=======
clear all; close all; addpath(genpath('./.'));
load('result-BSGS/GEN-Run_1_2017-05-07_14-37');
[Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
>>>>>>> origin/Standard
parm.k.covar = gen.covar;
parm.k.covar.range0 = fliplr(gen.covar.range0) ./ [grid_gen.dy grid_gen.dx];

<<<<<<< HEAD
parm.saveit = false;

parm.seed_path = 'shuffle';
parm.seed_search = 'shuffle';
parm.seed_U = 'shuffle';

parm.k.wradius = 3;
parm.k.lookup = false;
parm.k.nb = 30;
parm.mg=1;
parm.par_n=4;

=======
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
>>>>>>> origin/Standard

% use the log of hyd. cond.
hd = sampling_pt(struct('x',1:grid_gen.nx,'y',1:grid_gen.ny),log(K_true),1,4);
hd.d = Nscore.forward(hd.d);

f0=kern.prior ./ sum(kern.prior);
nx = grid_gen.nx;
ny = grid_gen.ny;




%% work on the weight
<<<<<<< HEAD
parm.aggr.method='cst';
parm.aggr.T = [0 .5 1]';

% parm.aggr.method='step';
% parm.aggr.T = [0 .5 1]';
=======
% parm.aggr.method='cst';
% parm.aggr.T = .5;

filename='ResStepVarious';
parm.aggr.method='step';
parm.aggr.T = [0  .001 .002 .003 .06 .1 .5 1]';
>>>>>>> origin/Standard

% parm.aggr.method='linear';
% parm.aggr.T = [ .1 .9];
    
% filename='ResSigmoid06';
% parm.aggr.method='sigmoid';
% parm.aggr.T = [ .06 Inf ; .06 1000; .06  100; .06  50; .06  20; .06  10];


parm.aggr.sum = 1;
<<<<<<< HEAD
parm.n_real  = parm.par_n*numel(parm.aggr.T)*2;
=======
parm.n_realisation  = parm.par_n*size(parm.aggr.T,1);
>>>>>>> origin/Standard


disp('Setup ok, Start')
Res =  BSS(nx,ny,hd,f0,sec,parm);
disp('Run finished')

<<<<<<< HEAD
id = grid_gen.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
Gamma_t_id = Gamma_t(id);
XY = kern.XY;
Sigma_d = Sigma.d(:);
dens = kern.dens(:);


E1 = nan(sum(id),parm.n_real);
E2 = nan(numel(kern.dens),parm.n_real);


for i_real=1:parm.n_real
   r = Res(:,:,i_real);
   gamma_x = variogram_gridded_perso(r);
   E1(:,i_real) = gamma_x(id)-Gamma_t_id;
   E2(:,i_real) = ksdensity([Sigma_d r(:)],XY)-dens;
=======
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
>>>>>>> origin/Standard
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
