%% Setting 
clear all
load('result-BSS/GEN-Run_1_2017-05-07_14-37');
[Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
addpath(genpath('./.'));
parm.k.covar = gen.covar;


parm.unit='';
parm.nscore = 1;
parm.par = 1;

parm.k.nb = [0 0 0 0 0; 10 10 10 10 20];
parm.cstk = 1;
parm.seed = 'shuffle';
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 1;
parm.plot.krig=0;
parm.path='quasirandom';
parm.path_random=1;
parm.kernel=kern;


parm.plot.krig = 0;
parm.par_n = 4;


Klog=K;
Klog.d=log(Klog.d);



% parm.n_realisation=0;parm.notify =0;
% [~,~,kern,~,~,~] =  BSGS(Klog,Sigma,grid_gen,parm);
% parm.notify=1;


Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
% [X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
% XY = [X(:),Y(:)];
%jpdf = ksdensity([Sigma.d(kern.id) Klog.d],XY);
%jpdf = ksdensity([Sigma.d(:) log(K_true(:))],XY);


%% work on the weight
% parm.aggr.method='AB';
% parm.aggr.a=0:.1:1;%.01:.01:.19;
% parm.aggr.b=5;
% [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
% parm.n_realisation  = parm.par_n*numel(parm.aggr.A)*4;
% figure(929);hold on;
% i_pt_temp=1:grid_gen.nxy;
% fx = @(a,b,x) (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
% for i_a=1:numel(parm.aggr.a)
%     a=parm.aggr.a(i_a);
%     for i_b=1:numel(parm.aggr.b)
%         b=parm.aggr.b(i_b);
%         plot(i_pt_temp,fx(a,b,i_pt_temp./grid_gen.nxy));
%     end
% end

parm.n_realisation  = parm.par_n*4;
parm.aggr.method='cst';
parm.aggr.w=1;

 
% parm.aggr.method='rad';
% parm.aggr.x=0:.05:1;%.01:.01:.19;
% parm.n_realisation  = parm.par_n*numel(parm.aggr.x)*4;
 

[Res,~,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);


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
    if strcmp(parm.aggr.method,'AB')
        parfor i_ab = 1:numel(parm.aggr.A)
            idAB=i_ab:numel(parm.aggr.A):parm.n_realisation;
            OF1(i_ab) = sqrt(mean(mean(E1(:,idAB),2).^2));
            OF2(i_ab) = sqrt(mean(mean(E2(:,idAB),2).^2));
        end
    elseif strcmp(parm.aggr.method,'rad')
        parfor i_x = 1:numel(parm.aggr.x)
            idx=i_x:numel(parm.aggr.x):parm.n_realisation;
            OF1(i_x) = sqrt(mean(mean(E1(:,idx),2).^2));
            OF2(i_x) = sqrt(mean(mean(E2(:,idx),2).^2));
        end
    else
        OF1 = sqrt(mean(mean(E1,2).^2));
        OF2 = sqrt(mean(mean(E2,2).^2));
    end
catch
    save('result-BSS/ResRAD_X_Norm.mat','Res','parm','kern','Nscore','E1','E2')
end

save('result-BSS/ResRAD_X_Norm.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')
%save('result-BSS/ResAB_A_1.mat','parm','E1','E2')
% save('C:\Users\Raphael\Desktop\to be compress\ResRAD_AB.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')


%% Calibration
% XY=[X(:),Y(:)];
% Sigma_d=Sigma.d(:);
% parm.n_realisation  = parm.par_n*1;
% parm.notify = 0;
% parm.aggr.fx = @(a,b,x) (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
% fun = @(x) OF_fx(Klog,Sigma,grid_gen,parm,id,Gamma_t,XY,jpdf,Sigma_d, x(1), x(2) );
% 
% 
% %Single min
% options = optimset('MaxFunEvals',3,'OutputFcn',@outfun);
% [x,fval,exitflag,output] = fminsearch(fun,[0 realmax],options);
% 
% save('result-BSS/Cal01.mat','parm','x','fval','exitflag','output');


% Many min
% options = optimset('Display','iter');
% [x,fval,exitflag,output] = fminsearch(fun,[.1 10],options)


% % Pareto Front
% nf = 2; % number of objective functions
% N = 5; % number of points for plotting
% x = zeros(N+1,1);
% f = zeros(N+1,nf);
% x0 = 0.5;
% goal=[0.01 .026];
% options = optimoptions('fgoalattain','Display','iter-detailed','MaxFunEvals',3);
% for r = 0:N
%     t = r/N; % 0 through 1
%     weight = [t,1-t];
%     if r~=0
%         x0=x(r,:);
%     end
%     [x(r+1,:),f(r+1,:),attainfactor,exitflag,output,lambda] = fgoalattain(fun,x0,goal,weight,[],[],[],[],[],[],[],options)
% end

