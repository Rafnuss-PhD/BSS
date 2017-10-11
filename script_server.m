
addpath(genpath('./.'));
load('result-BSS/GEN-Run_1_2017-05-07_14-37');

Nscore = nscore(kern, struct('nscore', 1), 0); %Prim.d, kern.axis_prim, 'pchip', 'pchip', parm.plot.ns

[~,s_id]=min(bsxfun(@minus,kern.axis_sec,Sigma.d(:)).^2,[],2);
sec_pdf = kern.dens(:,s_id);
sec.pdf = bsxfun(@times, sec_pdf, 1./sum(sec_pdf));
sec.axis = Nscore.forward(kern.axis_prim);

parm.k.covar = gen.covar;
parm.k.covar.range0 = fliplr(gen.covar.range0) ./ [grid_gen.dy grid_gen.dx];


parm.saveit = false;

parm.seed_path = 'shuffle';
parm.seed_search = 'shuffle';
parm.seed_U = 'shuffle';

parm.k.wradius = 3;
parm.k.lookup = false;
parm.k.nb = 30;
parm.mg=1;


% use the log of hyd. cond.
hd = sampling_pt(struct('x',1:grid_gen.nx,'y',1:grid_gen.ny),log(K_true),1,4);
hd.d = Nscore.forward(hd.d);

f0=kern.prior ./ sum(kern.prior);
nx = grid_gen.nx;
ny = grid_gen.ny;




%% work on the weight

parm.aggr.method='cst';
parm.aggr.T = 0;

% parm.aggr.method='step';
% parm.aggr.T = [0 .5 1]';

% parm.aggr.method='cst';
% parm.aggr.T = .5;



% parm.aggr.method='linear';
% parm.aggr.T = [ .1 .9];
    
% parm.aggr.method='sigmoid';
% parm.aggr.T = [ .06 Inf ; .06 1000; .06  100; .06  50; .06  20; .06  10];


filename='ResCst0';
parm.aggr.sum = 1;

parm.par_n=48;
parpool(parm.par_n)
parm.n_real  = parm.par_n*numel(parm.aggr.T);


disp('Setup ok, Start')
n=6;
Res=nan(ny,nx,parm.n_real*n);
for i=1:n
    ii = (i-1)*parm.n_real + (1:parm.n_real);
    Res(:,:, ii ) =  BSS(nx,ny,hd,f0,sec,parm);
    disp(['Sim:' num2str(i/n)])
end
disp('Run finished')



id = grid_gen.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
Gamma_t_id = Gamma_t(id);
XY = kern.XY;
Sigma_d = Sigma.d(:);
dens = kern.dens(:);

E1 = nan(sum(id),parm.n_real);
E2 = nan(numel(kern.dens),parm.n_real);


parfor i_real=1:parm.n_real
   r = Res(:,:,i_real);
   gamma_x = variogram_gridded_perso(r);
   E1(:,i_real) = gamma_x(id)-Gamma_t_id;
   E2(:,i_real) = ksdensity([Sigma_d r(:)],XY)-dens;
end
disp('Error Computed')


try
    clear OF1 OF2
    parfor i_t = 1:size(parm.aggr.T,1)
        ii_t=i_t:size(parm.aggr.T,1):parm.n_real;
        OF1(i_t) = sqrt(mean(mean(E1(:,ii_t),2).^2));
        OF2(i_t) = sqrt(mean(mean(E2(:,ii_t),2).^2));
    end
catch
    %save(['result-BSGS/' filename],'Res','parm','kern','Nscore','E1','E2')
end

save(['result-BSGS/' filename],'parm','E1','E2','OF1','OF2')
disp('file written')
