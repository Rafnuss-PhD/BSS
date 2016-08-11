
%% Create a standard case

% Grid size
gen.xmax = 240; %total length in unit [m]
gen.ymax = 20; %total hight in unit [m]

% Scale define the subdivision of the grid (multigrid). At each scale, the
% grid size is $(2^gen.sx+1) \times (2^gen.sy+1)$ 
gen.sx = 10;
gen.sy = 7;

% Generation Method: All method generate with FFTMA a gaussian field.
% 'Normal'              with normal distribution \sim N(gen.mu,gen.std)
% 'LogNormal'   
% 'fromRho':            log transform it with the parameter defined below 
% 'fromK':              generate with FFTMA a field of Hyraulic conductivity and log transform it with the parameter defined below 
gen.method              = 'fromRho';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar.modele        = [4 70 7 0];%; 1 1 1 0]; % covariance structure
gen.covar.c             = [1];%; 0.01]; 
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 256;
gen.Rho.grid.ny           = 32; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'Calibration-';
gen.seed                = 123456;

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);


%% Create simulation

clear all;clc; dbstop if error
familyname='4-Calibration/';

load(['result-time-log-calibration/' familyname 'GEN-Calibration-_2016-06-15_18-14']);
parm.gen=gen;


%parm.w_krig = @(x,a,b)  (atan(a*b) - atan(b*(a - x)))/(atan(a*b) - atan(b*(a - 1)));
a=.5;
b=100;


parm.w_krig = @(x) (atan(a*b) - atan(b*(a - x)))/(atan(a*b) - atan(b*(a - 1)));
%hold on; plot(0:.01:1,parm.w_krig(0:.01:1));
%parm.w_sec = @(x) .4;

% parm.scale=[1:8;1:7 7];

parm.n_realisation = 40;
parm.par = 1;

[~,~,~,~,~, filename] = BSGS(sigma,Sigma,sigma_true,grid_gen,parm);
% filename = 'SIM-_2016-03-15_16-37-50';




%%


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


