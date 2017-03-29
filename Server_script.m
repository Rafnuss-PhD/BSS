disp('start')
%% Generation

% Grid size
gen.xmax = 200; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.sx = 9;
gen.sy = 6;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'fromPhi';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar(1).model      = 'spherical';
gen.covar(1).range0     = [10 3];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 200;
gen.Rho.grid.ny           = 15; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 2; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'test';
gen.seed                = 'default';

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);


%% Setting 


% load('result-BSS/GEN-test_2017-01-11_10-59');
% addpath(genpath('./.'));
% parm.k.covar = gen.covar;
% 
% 
% parm.unit='';
% parm.nscore = 1;
% parm.par = 1;
% 
% parm.k.nb = [0 0 0 0 0; 10 10 10 10 20];
% parm.cstk = 1;
% parm.seed = 'shuffle';
% parm.scale=[1:grid_gen.sx;1:grid_gen.sy 5 5 5]; % no multigrid
% parm.saveit = false;
% parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
% parm.k.quad = 0;
% parm.k.wradius = 1.3;
% parm.plot.kernel=0;
% parm.plot.ns= 0;
% parm.plot.krig=0;
% parm.path='quasirandom';
% parm.path_random=1;
% 
% parm.plot.krig=0;
% parm.par_n =44;
% 
% 
% Klog=K;
% Klog.d=log(Klog.d);
% 
% parm.n_realisation=0;parm.notify =0;
% [~,t,kern,~,~,~] =  BSGS(Klog,Sigma,grid_gen,parm);
% parm.notify=1;
% 
% id = grid_gen.x<parm.k.covar(1).range(1)*parm.k.wradius;
% Gamma_t = (1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)))';
% [X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
% jpdf = ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]);


%% work on the weight
% parm.aggr.a=[0:.1:.4];
% parm.aggr.b=[50];
% parm.aggr.fx = @(a,b,x) (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
% [parm.aggr.A, parm.aggr.B] = meshgrid(parm.aggr.a, parm.aggr.b);
% parm.aggr.w=0.5;
% 
% parm.n_realisation  = parm.par_n*numel(parm.aggr.A)*6;
% 
% % figure(99);hold on;
% % i_pt_temp=1:grid_gen.nxy;
% % for i_a=1:numel(parm.aggr.a)
% %     a=parm.aggr.a(i_a);
% %     for i_b=1:numel(parm.aggr.b)
% %         b=parm.aggr.b(i_b);
% %         plot(i_pt_temp,parm.aggr.fx(a,b,i_pt_temp./grid_gen.nxy));
% %     end
% % end
% [Res,t,kern,~,~,Nscore] =  BSGS(Klog,Sigma,grid_gen,parm);
% 
% 
% 
% E1 = nan(sum(id),parm.n_realisation);
% E2 = nan(numel(jpdf),parm.n_realisation);
% 
% parfor i_realisation=1:parm.n_realisation
%     gamma_x = variogram_gridded_perso(Res.m_ns{i_realisation});
%     E1(:,i_realisation) = gamma_x(id)-Gamma_t;
%     E2(:,i_realisation) = ksdensity([Sigma.d(:) Res.m{i_realisation}(:)],[X(:),Y(:)])-jpdf;
% end
% 
% parfor i_ab = 1:numel(parm.aggr.A)
%     idAB=i_ab:numel(parm.aggr.A):parm.n_realisation;
%     OF1(i_ab) = sqrt(mean(mean(E1(:,idAB),2).^2));
%     OF2(i_ab) = sqrt(mean(mean(E2(:,idAB),2).^2));
% end
% 
% save('result-BSS/ResAB5.mat','Res','parm','kern','Nscore','E1','E2','OF1','OF2')


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

