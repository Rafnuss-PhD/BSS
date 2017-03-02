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
parm.k.nb = [0 0 0 0 0; 10 10 10 10 100];
parm.cstk = true;
parm.seed = 'shuffle';
parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
parm.saveit = false;
parm.k.method = 'sort'; % sort, sbss (29.6) or smart minkmex
parm.k.quad = 0;
parm.k.wradius = 1.3;
parm.plot.kernel=0;
parm.plot.ns= 0;
parm.path='quasirandom';
parm.path_random='false';

parm.plot.krig=0;
Klog=K;
Klog.d=log(Klog.d);

% work on the weight
a=.1;
b=10;
parm.aggr.fx = @(parm,grid,i_scale,i_pt) 0;%(atan(a*b) - atan(b*(a -  i_pt./grid{end}.nxy )))/(atan(a*b) - atan(b*(a - 1)));

% figure(99);hold on; 
% i_pt_temp=1:grid_gen.nxy;
% plot(i_pt_temp,parm.aggr.fx(parm,grid_gen,1,i_pt_temp));

[Res, t, kern, k, ~, Nscore, filename] =  BSGS(Klog,Sigma,grid_gen,parm);





%%
% Histogram:
figure(1);clf; hold on;
for i_realisation=1:parm.n_realisation
    ksdensity(Res{end}.m_ns{i_realisation}(:))
    [~,p(i_realisation)] = adtest(Res{end}.m_ns{i_realisation}(:),'Distribution',makedist('normal','mu',0,'sigma',1));
end
a=-3:.01:3;
plot(a,normpdf(a,0,1),'--k','LineWidth',2)
legend(strread(num2str(p),'%s'))



% Variogram
for i_realisation=1:parm.n_realisation
    [gamma_x_y{i_realisation}, gamma_y_y{i_realisation}] = variogram_gridded_perso(Res{end}.m_ns{i_realisation});
end

figure(2);clf; 
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
id_prim = std(Res2D).^2<eps; % remove hard data
Cc = C;
Cc(:,id_prim)=[];
Cc(id_prim,:)=[];
Res2Dc=Res2D;
Res2Dc(:,id_prim)=[];
Klogkrig2 = Klogkrig(:);
Klogkrig2(id_prim)=[];
ResM=mean(Res2Dc);
ResC=cov(Res2Dc);

T2 = (ResM-Klogkrig2(:)')*inv(Cc/parm.n_realisation)*(ResM-Klogkrig2(:)')';

P=1-chi2cdf(T2,size(Res2Dc,2));

figure(4); clf;
subplot(1,3,1); imagesc(ResC); colorbar;
subplot(1,3,2); imagesc(Cc); colorbar;
subplot(1,3,3); imagesc(ResC-Cc); colorbar;

a=reshape(sum(cov(Res2D)-C),grid_gen.ny, grid_gen.nx);
figure; imagesc(a)

%% OLD

OF_true.mean = mean(Prim_true(:));
OF_true.std = std(Prim_true(:));
% Variogram
[OF_true.vario_x_field, OF_true.vario_y_field] = variogram_gridded_perso(Prim_true);
OF_true.vario_x_theo = semivariogram1D(grid{end}.x,1,parm.gen.covar.modele(1,2),parm.gen.covar.modele(1),0);
OF_true.vario_y_theo = semivariogram1D(grid{end}.y,1,parm.gen.covar.modele(1,3),parm.gen.covar.modele(1),0);


of(1)=OF('SIM-00_2016-06-16_16-44-10');
of(2)=OF('SIM-025_2016-06-17_09-43-54');
of(3)=OF('SIM-05_2016-06-17_12-34-45');
of(4)=OF('SIM-10_2016-06-16_16-12-32');


name={'0','.25','.5','1'};
col={'b','r','g','m'};


figure(101)
subplot(2,2,[1 2]); hold on;
for i=1:numel(of)
    plot(of(i).hist_pts, mean(of(i).hist), col{i})
    plot(of(i).hist_pts, mean(of(i).hist)-std(of(i).hist) ,['--' col{i}])
    plot(of(i).hist_pts, mean(of(i).hist)+std(of(i).hist) ,['--' col{i}])
end
plot(-4:.01:4,normpdf(-4:.01:4,0,1),'--k')
axis tight;
legend(name)

subplot(2,2,3); hold on;
plot([1 numel(of)],[0 0],'--k')
boxplot(reshape([of.mean],40,numel(of)),name)

subplot(2,2,4); hold on;
plot([1 numel(of)],[1 1],'--k')
boxplot(reshape([of.std],40,numel(of)),name)



figure(102)
subplot(1,2,1); hold on;
plot(grid{end}.x, OF_true.vario_x_theo,'--k')
plot(grid{end}.x, OF_true.vario_x_field,'--k')
for i=1:numel(of)
    plot(grid{end}.x, mean(of(i).gamma_x,2),'.-')
end
axis tight; xlim([0 parm.gen.covar.modele(1,2)*1.3])
legend({'theory','true field',name{:}})
subplot(1,2,2); hold on;
plot(grid{end}.y, OF_true.vario_y_theo,'--k')
plot(grid{end}.y, OF_true.vario_y_field,'--k')
for i=1:numel(of)
    plot(grid{end}.y, mean(of(i).gamma_y,2),'-')
end
axis tight; xlim([0 parm.gen.covar.modele(1,3)*1.3])
legend({'theory','true field',name{:}})


figure(103);
c_axis=[min(parm.gen.Rho.f.output.pseudo_interp(:)), max(parm.gen.Rho.f.output.pseudo_interp(:))];
subplot(numel(of)+1,1,1);
imagesc(parm.gen.Rho.f.output.pseudo_interp); caxis(c_axis); axis tight;
for i=1:numel(of)
    subplot(numel(of)+1,1,i+1);
    imagesc(mean(of(i).ert_interp,3)); caxis(c_axis); axis tight;
    xlabel(name{i})
end

figure(104); hold on
for i=1:numel(of)
    plot(of(i).w(1,:),[col{i} 'x-'])
end


figure(105); hold on;
plot([1 numel(of)],[0 0],'--k')
boxplot(reshape([of.ert_err],40,numel(of)),name)

figure(106); hold on;
plot([1 numel(of)],[0 0],'--k')
boxplot(reshape([of.H],40,numel(of)),name)


