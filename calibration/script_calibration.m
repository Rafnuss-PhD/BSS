clc; % clear the command line
addpath(genpath('./.'));  % Add folder and sub-folder to path
dbstop if error  % activate debug in error

%% Create a standard case

% Grid size
gen.xmax = 200; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.sx = 8;
gen.sy = 5;

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
gen.covar(1).range0     = [15 4];
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
gen.Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'Calibration-';
gen.seed                = 'default';

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);


%% Create simulation

clear all; clc;  close all; dbstop if error
load('result-BSS/GEN-test_2017-01-11_10-59');


parm.k.covar = gen.covar;

parm.unit='';
parm.nscore = 1;
parm.par = 0;
parm.n_realisation  = 100;
parm.k.nb = [0 0 0 0 0; 10 10 10 10 20];
parm.cstk = true;
parm.seed = 'default';
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 0;
parm.plot.krig=0;
parm.path='quasirandom';
parm.path_random=0;
parm.notify=1;

Klog=K;
Klog.d=log(Klog.d);

% work on the weight
a=.1;
b=10;

parfor i=1:11
    parm1=parm;
    parm1.aggr.fx = @(parm,grid,i_scale,i_pt) (i-1)/10 ;%(atan(a*b) - atan(b*(a -  i_pt./grid{end}.nxy )))/(atan(a*b) - atan(b*(a - 1)));
    
    % figure(99);hold on;
    % i_pt_temp=1:grid_gen.nxy;
    % plot(i_pt_temp,parm.aggr.fx(parm,grid_gen,1,i_pt_temp));
    
    Res{i} =  BSGS(Klog,Sigma,grid_gen,parm1);

end




%% TEST HYPOTHESE
% Histogram:
figure(1);clf; hold on;
for i_realisation=1:parm.n_realisation
    ksdensity(Res1{end}.m_ns{i_realisation}(:))
    [~,p(i_realisation)] = adtest(Res1{end}.m_ns{i_realisation}(:),'Distribution',makedist('normal','mu',0,'sigma',1));
end
a=-3:.01:3;
plot(a,normpdf(a,0,1),'--k','LineWidth',2)
legend(strread(num2str(p),'%s'))

clear gamma_x

% Variogram
id = grid_gen.x<parm.k.covar(1).range(1);
Gamma_t = [(1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)))];%...
    %(1-parm.k.covar(1).g(grid_gen.y(id2)/parm.k.covar(1).range(2)))];
gamma_x=nan(grid_gen.nx,parm.n_realisation);
parfor i_realisation=1:parm.n_realisation
    gamma_x(:,i_realisation) = variogram_gridded_perso(Res1{end}.m_ns{i_realisation});
end

err= bsxfun(@min,gamma_x(id,:),Gamma_t');
figure; hold on;
for i_realisation= 1:300
    E1(i_realisation) = sqrt(mean((mean(err(:,randi(parm.n_realisation,i_realisation)),2)).^2));
end
plot(E1)

X = [gamma_x(id,:)'];% gamma_y(id2,:)'];
[n,p]=size(X);
Gamma = mean(X);
S_Gamma = cov(X,1);

T2 = (n-1) * (Gamma-Gamma_t) / S_Gamma * (Gamma-Gamma_t)';
% alpha=.05;
%T2max = p*(n-1)/(n-p)*finv(1-alpha,n,n-p);
P=1-fcdf((n-p)/((n-1)*p)*T2, p, n-p);  %Probability that null Ho: is true. 5.4252e-05

figure; hold on;
plot(X','-k')
plot(Gamma,'-k','linewidth',5)
plot(Gamma_t,'--r','linewidth',2)
legend(num2str(T2),num2str(P))



figure;clf; 
subplot(1,2,1); hold on
id= grid_gen.x<parm.k.covar(1).range(1)*1.2;
for i_realisation=1:parm.n_realisation
    h1=plot(grid_gen.x(id),gamma_x_y{i_realisation}(id),'Color', [.4 .4 .4]);
end
h2=plot(grid_gen.x(id),(1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1))),'--k','linewidth',3);
legend([h2, h1],'Theorical model C(h)','simulation(s)','Location','northwest')
ylabel('Horizontal'); xlabel('m');axis tight

subplot(1,2,2); hold on
id= grid_gen.y<parm.k.covar(1).range(2)*1.5;
for i_realisation=1:parm.n_realisation
    h1=plot(grid_gen.y(id),gamma_y_y{i_realisation}(id),'Color', [.4 .4 .4]);
end
h2=plot(grid_gen.y(id),(1-parm.k.covar(1).g(grid_gen.y(id)/parm.k.covar(1).range(2))),'--k','linewidth',3);
legend([h2, h1],'Theorical model C(h)','simulation(s)','Location','northwest')
ylabel('Vertical'); xlabel('m'); axis tight


% Expected field
a0_C=covardm_perso([Klog.x Klog.y],[grid_gen.X(:) grid_gen.Y(:)],parm.k.covar);
ab_C=covardm_perso([Klog.x Klog.y],[Klog.x Klog.y],parm.k.covar);
Klogkrig=reshape((ab_C \ a0_C)'*Nscore.forward(Klog.d),grid_gen.ny,grid_gen.nx);

Res3Dmean = mean(Res3D,3);
figure(3); clf;
subplot(3,1,1); imagesc(Klogkrig)
subplot(3,1,2); imagesc(Res3Dmean)
subplot(3,1,3); imagesc(Res3Dmean-Klogkrig)


% Covariances
C=sparse(covardm_perso([grid_gen.X(:) grid_gen.Y(:)],[grid_gen.X(:) grid_gen.Y(:)],parm.k.covar));
Res3D = reshape([Res{end}.m_ns{:}], [grid_gen.ny, grid_gen.nx, parm.n_realisation]);
Res2D=reshape(Res3D,[grid_gen.nxy, parm.n_realisation])';

ResC = cov(Res2D);

figure; clf;
% subplot(1,3,1); imagesc(ResC); colorbar;
% subplot(1,3,2); imagesc(C); colorbar;
% subplot(1,3,3); 
imagesc(ResC-C); colorbar;


% Secondary variable
figure; clf;
subplot(3,3,[1 3]); imagesc(Sigma.d); colorbar;
subplot(3,3,[4 6]); imagesc(Res3Dmean); colorbar;

[X,Y] = meshgrid(linspace(min(Sigma.d(:)),max(Sigma.d(:)),100), linspace(min(Res{1}.m{1}(:)),max(Res{1}.m{1}(:)),100));
j1 = ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]);
j1 = reshape(j1,size(X,1),size(X,2));
j2=nan(size(X,1),size(X,2), parm.n_realisation);
for i_realisation=1:parm.n_realisation
    j2(:,:,i_realisation) = reshape(ksdensity([Sigma.d(:) Res{1}.m{i_realisation}(:)],[X(:),Y(:)]),size(X,1),size(X,2));
end

subplot(3,3,7); imagesc(j1)
subplot(3,3,8); imagesc(mean(j2,3))
subplot(3,3,9); imagesc(mean(j2,3)-j1)



j1=[sigma.d(:) Klog.d];
j2=[];
for i_realisation=1:1%parm.n_realisation
    j2 = [j2; Sigma.d(:) Res{1}.m{i_realisation}(:)];
end

[p,D] = kstest2d(j1,j2);


figure; hold on;
subplot(1,2,1); ksdensity(j1)
subplot(1,2,2); ksdensity(j2)

%% AB
clear E1 E2
id = grid_gen.x<parm.k.covar(1).range(1)*parm.k.wradius;
Gamma_t = [(1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)))];%...
    %(1-parm.k.covar(1).g(grid_gen.y(id2)/parm.k.covar(1).range(2)))];
gamma_x=nan(grid_gen.nx,parm.n_realisation);
parfor i_realisation=1:parm.n_realisation
    gamma_x(:,i_realisation) = variogram_gridded_perso(Res{end}.m_ns{i_realisation});
end
err= gamma_x(id,:)-repmat(Gamma_t',1,size(gamma_x,2));
E1 = sqrt(mean(err.^2));
for i_ab = 1:numel(parm.aggr.A)
    idAB(i_ab)=i_ab:numel(parm.aggr.A):parm.n_realisation;
    E1AB(i_ab) = sqrt(mean(mean(err(:,idAB(i_ab)),2).^2));
end


[X,Y] = meshgrid(linspace(min(sigma.d(:)),max(sigma.d(:)),100),...
    linspace(min(Klog.d(:)),max(Klog.d(:)),100));
%jpdf = reshape(ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]),size(X,1),size(X,2));
%Res3D = reshape([Res{end}.m{:}], [grid_gen.ny, grid_gen.nx, parm.n_realisation]);
jpdf = ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]);
Res3D = [Res{end}.m{:}];
parfor i_realisation=1:parm.n_realisation
    a = Res3D(:,:,i_realisation);
    Resjpdf=reshape(ksdensity([Sigma.d(:) a(:)],[X(:),Y(:)]),size(X,1),size(X,2));
    E2(i_realisation) = sqrt(mean((Resjpdf(:) - jpdf(:)).^2));
end

for i_ab = 1:numel(parm.aggr.A)
    idAB{i_ab}=i_ab:numel(parm.aggr.A):parm.n_realisation;
    a = [Res{end}.m{idAB{i_ab}}];
    Resjpdf=ksdensity([Sigma.d(:) a(:)],[X(:),Y(:)]);
    E2(i_realisation) = sqrt(mean((Resjpdf(:) - jpdf(:)).^2));
end



figure;
subplot(2,1,1); plot(OF1);
subplot(2,1,2); plot(OF2);

figure; hold on;
for i=1:size(parm.aggr.A,1)
    for j=1:size(parm.aggr.A,2)
        id = sub2ind(size(parm.aggr.A),i,j);
        ii = id:parm.par_n:parm.n_realisation;
        plot(E1(ii),E2(ii),'o')
        keyboard
    end
end


figure; hold on;
for i=1:numel(parm.aggr.w)
    ii = i:numel(parm.aggr.w):parm.n_realisation;
    plot(E1(ii),E2(ii),'o')
end

figure; hold on;
for i=1:numel(parm.aggr.w)
    ii = i:numel(parm.aggr.w):parm.n_realisation;
    subplot(1,2,1); hold on; histogram(E1(ii))
    subplot(1,2,2); hold on; histogram(E2(ii))
end


%% Error
id = grid_gen.x<parm.k.covar(1).range(1);
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)));
[X,Y] = meshgrid(linspace(min(sigma.d(:)),max(sigma.d(:)),100),...
    linspace(min(Klog.d(:)),max(Klog.d(:)),100));
jpdf = reshape(ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]),size(X,1),size(X,2));
Sigmad3D = repmat(Sigma.d,1,1,parm.n_realisation);

clear gamma_x Resjpdf E1 E2 Gamma Resjpdf
parfor i=1:numel(Res)
    % Variogram
    gamma_x=nan(grid_gen.nx,parm.n_realisation);
    for i_realisation=1:parm.n_realisation
        [gamma_x(:,i_realisation), ~] = variogram_gridded_perso(Res{i}{end}.m_ns{i_realisation});
    end
    Gamma{i} = mean(gamma_x(id,:),2)';
    E1(i) = sqrt(mean((Gamma{i} - Gamma_t).^2));

    % Joint PDF
    Res3D = reshape([Res{i}{end}.m{:}], [grid_gen.ny, grid_gen.nx, parm.n_realisation]);
    Resjpdf{i}=reshape(ksdensity([Sigmad3D(:) Res3D(:)],[X(:),Y(:)]),size(X,1),size(X,2));
    E2(i) = sqrt(mean((Resjpdf{i}(:) - jpdf(:)).^2));
end


figure;
subplot(1,2,1); plot(E1,'x'); xlabel('Variogram Error')
subplot(1,2,2); hold on; plot(E2,'x'); xlabel('Joint PDF')

figure;
scatter(E1,E2)

figure; hold on;
for i=1:numel(Res)
    plot(Gamma{i});
end
plot(Gamma_t,'-r','linewidth',3)


%% Calibration
parm.n_realisation  = 30;
parm.par = 1;
parm.seed = 'default';

id = grid_gen.x<parm.k.covar(1).range(1);
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)));
[X,Y] = meshgrid(linspace(min(sigma.d(:)),max(sigma.d(:)),100),...
    linspace(min(Klog.d(:)),max(Klog.d(:)),100));
jpdf = reshape(ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]),size(X,1),size(X,2));
Sigmad3D = repmat(Sigma.d,1,1,parm.n_realisation);


fun = @(x) OF(Klog,Sigma,grid_gen,parm,id,Gamma_t,X,Y,jpdf,Sigmad3D, x);
fun(0)
fun(1)


% Pareto Front
nf = 2; % number of objective functions
N = 5; % number of points for plotting
onen = 1/N;
x = zeros(N+1,1);
f = zeros(N+1,nf);
x0 = 0.5;
goal=[0.01 .028]; % 0.942573   0.00879023 
options = optimoptions('fgoalattain','Display','iter','PlotFcn',{@optimplotx @optimplotfval});
for r = 0:N
    t = onen*r; % 0 through 1
    weight = [t,1-t];
    [x(r+1,:),f(r+1,:)] = fgoalattain(fun,x0,goal,weight,...
        [],[],[],[],[],[],[],options);
end

figure
plot(f(:,1),f(:,2),'k.');
xlabel('f_1')
ylabel('f_2')



x0= .5*ones(grid.)

[xf,fval,exitflag,output,grad,hessian] = fminbnc( @(x) OF(Klog,Sigma,sigma,grid_gen,parm,x), x0,options)



