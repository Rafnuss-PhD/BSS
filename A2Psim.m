%% SCRIPT_BSS.m
% This file contain the script which can run the entire algorithm. It is
% where the main options/parameters are set and the main fonctions are
% lauch.
%
% # DATA CREATION : Data are created either by reading mat file or
% generated with physical relationship. 
% # BSS : Baysian Sequential Simulation
% # PLOT : Graphical visualisation of data
%
%
% Variable naming convension:
% * rho     : Electrical resisitivity [\omega.m]
% * sigma   : Electrical conductivity = 1/rho [mS/m]
% * phi     : Porosity [-]
% * K       : Hydraulic conductivity [m/s]
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)


clc; % clear the command line
addpath(genpath('./.'));  % Add folder and sub-folder to path
dbstop if error  % activate debug in error

%% DATA CREATION
% This section gather all possible way to create the data. |gen| struct
% store the parameter and |data_generation.m| compute everything

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
gen.method              = 'fromPhi';

% Generation parameter
gen.samp                = 1;          % Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
gen.samp_n              = 4;          % number of well or number of point
gen.covar(1).model      = 'exponential';
gen.covar(1).range0     = [27 2.7];
gen.covar(1).azimuth    = 0;
gen.covar(1).c0         = 1;
gen.covar               = kriginginitiaite(gen.covar);
gen.mu                  = 0.27; % parameter of the first field. 
gen.std                 = .05;
gen.Rho.method          = 'R2'; % 'Paolo' (default for gen.method Paolo), 'noise', 'RESINV3D'

% Electrical inversion
gen.Rho.grid.nx           = 240;
gen.Rho.grid.ny           = 20; % log-spaced grid.
gen.Rho.elec.spacing      = 2; % in grid spacing unit.
gen.Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
gen.Rho.dmin.res_matrix   = 2; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
gen.Rho.dmin.tolerance    = 10;

% Other parameter
gen.plotit              = true;      % display graphic or not (you can still display later with |script_plot.m|)
gen.saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
gen.name                = 'test';
gen.seed                = 'default';

% Run the function
data_generation(gen);
%[fieldname, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen);

%% Test Areal-to-point Kriging simulation
covar.model='exponential'; covar.range0=[20 20]; covar.c0=1;covar.azimuth=0;
gen.covar= kriginginitiaite(covar);
grid_gen.x=1:100;
grid_gen.y=1:100;
[grid_gen.X,grid_gen.Y] = meshgrid(grid_gen.x, grid_gen.y);
grid_gen.nxy =numel(grid_gen.X);
grid_gen.nx=numel(grid_gen.x);
grid_gen.ny=numel(grid_gen.y);
zt = fftma_perso(gen.covar, grid_gen);
NSigma.dx=5;
NSigma.dy=5;
NSigma.x = NSigma.dx/2:NSigma.dx:grid_gen.x(end);
NSigma.y = NSigma.dy/2:NSigma.dy:grid_gen.y(end);
[NSigma.X,NSigma.Y] = meshgrid(NSigma.x, NSigma.y);

% Built G
G = zeros(numel(NSigma.X),grid_gen.nxy);
for ij=1:grid_gen.nxy
     [~,id_z]=min((grid_gen.X(ij)-NSigma.X(:)).^2 + (grid_gen.Y(ij)-NSigma.Y(:)).^2);
     G(id_z,ij)=1;%/(NSigma.dx*NSigma.dy);
end
for ij = 1:numel(NSigma.X)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end

% Compute d
d = reshape(G * zt(:), numel(NSigma.y), numel(NSigma.x));

% kriging prediction
Cz = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [grid_gen.X(:) grid_gen.Y(:)], gen.covar);
Czd = Cz * G';
Cd = G * Czd;
W=zeros(grid_gen.nxy,numel(NSigma.X));
for ij=1:grid_gen.nxy
    W(ij,:) = Cd \ Czd(ij,:)';
end
zh = reshape( W * d(:), grid_gen.ny,grid_gen.nx);

% Build unconditional realization
zs = fftma_perso(gen.covar, grid_gen);
ds = reshape(G * zs(:), numel(NSigma.y), numel(NSigma.x));
zhs = reshape( W * ds(:), grid_gen.ny,grid_gen.nx);

% Build condtional realization
zcs = zh + (zs - zhs);
dcs = reshape(G * zcs(:), numel(NSigma.y), numel(NSigma.x));


% Figure
figure(55);clf; 
c_axis=[ min(zt(:)) max(zt(:))];
subplot(3,3,1); imagesc(grid_gen.x, grid_gen.y, zt); caxis(c_axis); title('zt True field');
subplot(3,3,2); imagesc(grid_gen.x, grid_gen.y, d); caxis(c_axis); title('d True field');
subplot(3,3,3); imagesc(grid_gen.x, grid_gen.y, zh); caxis(c_axis); title('zh Kriging prediction based on the known averaged');
subplot(3,3,4); imagesc(grid_gen.x, grid_gen.y, zs); caxis(c_axis); title('zs Unconditional Realization');
subplot(3,3,5); imagesc(grid_gen.x, grid_gen.y, ds); caxis(c_axis); title('ds Unconditional Realization');
subplot(3,3,6); imagesc(grid_gen.x, grid_gen.y, zhs); caxis(c_axis); title('zhs Kriging prediction based on the unconditional realization');
subplot(3,3,7); imagesc(grid_gen.x, grid_gen.y, zcs); caxis(c_axis); title('zcs Conditionale Realization');
subplot(3,3,8); imagesc(grid_gen.x, grid_gen.y, dcs); caxis(c_axis); title('dcs Conditional Realization');
subplot(3,3,9); imagesc(grid_gen.x, grid_gen.y, dcs-d);  title('Error'); caxis(c_axis);


figure(66); clf
id_z = grid_gen.x < gen.covar.range0(1)*1.5;
[gamma_zt]=variogram_gridded_perso(zt);
id_d = NSigma.x < gen.covar.range0(1)*1.5;
[gamma_d]=variogram_gridded_perso(d);
subplot(2,2,1); hold on; plot(grid_gen.x(id_z)-grid_gen.x(1), 1-Cz(1,id_z)); plot(grid_gen.x(id_z)-grid_gen.x(1), gamma_zt(id_z));
subplot(2,2,2); hold on; plot(NSigma.x(id_d)-NSigma.x(1), Cd(1)- Cd(1,id_d)); plot(NSigma.x(id_d)-NSigma.x(1), gamma_d(id_d));
subplot(2,2,3); histogram(zt(:));
subplot(2,2,4); histogram(d(:));


%% Areal-to-point kriging simulation with Actual Data
clear all; close all
load('result-A2PK/GEN-Run_1_2017-05-07_14-37');


% Nscore of Sigma. 
[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
NSigma.x=Sigma.x_raw; NSigma.y=Sigma.y_raw; 
[NSigma.X, NSigma.Y] = meshgrid(Sigma.x_raw, Sigma.y_raw);
d = reshape(Nscore.forward(Sigma.d_raw(:)) ,numel(NSigma.y),numel(NSigma.x));
zt = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);

% Built G
G = zeros(numel(NSigma.X),grid_gen.nxy);
for ij=1:grid_gen.nxy
     [~,id_z]=min((grid_gen.X(ij)-NSigma.X(:)).^2 + (grid_gen.Y(ij)-NSigma.Y(:)).^2);
     G(id_z,ij)=1;%/(NSigma.dx*NSigma.dy);
end
% % Different version to use Sigma.rad
% for ij = 1:numel(NSigma.X)
%     G(ij,G(ij,:)==1) = Sigma.rad_raw(ij) ./ sum(G(ij,:)==1);
%     G(ij,G(ij,:)==0) = ( 1-Sigma.rad_raw(ij)) ./ sum(G(ij,:)==0);
% end
for ij = 1:numel(NSigma.X)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end


c_axis=[ min(zt(:)) max(zt(:)) ];
figure(22); clf;
subplot(2,1,1); imagesc(grid_gen.x, grid_gen.y, zt); caxis(c_axis); title('True field')
subplot(2,1,2); surf(NSigma.x, NSigma.y, d,'EdgeColor','none'); caxis(c_axis); title('Inverted field'); view(2); axis tight; set(gca,'Ydir','reverse'); box on


% Compare actual ERT inverted and corresponding averaged Electricit
dG = reshape(G * zt(:), numel(NSigma.y), numel(NSigma.x));
Gt = d(:) * zt(:)' / (zt(:) * zt(:)');
dGt = reshape(Gt * zt(:), numel(NSigma.y), numel(NSigma.x));

figure(44);
subplot(3,1,1); imagesc(NSigma.x, NSigma.y, d); caxis(c_axis); title('Inverted field')
subplot(3,1,2); imagesc(NSigma.x, NSigma.y, dG); caxis(c_axis); title('Averaging true field G * z_{true}')
subplot(3,1,3); imagesc(NSigma.x, NSigma.y, dGt); caxis(c_axis); title('Averaging true field Gt * z_{true}')

% kriging prediction
Cz = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [grid_gen.X(:) grid_gen.Y(:)], gen.covar);
tic;
Czd = Cz * G';
toc; tic
Cd = G * Czd;
toc;
W=zeros(grid_gen.nxy,numel(NSigma.X));
parfor ij=1:grid_gen.nxy
    W(ij,:) = Cd \ Czd(ij,:)';
end
zh = reshape( W * d(:), grid_gen.ny,grid_gen.nx);

 % save('result-A2PK/A2Psim-Run_1_2017-05-07_14-37_2','Cd','Czd','W','G','Gt','d','zt','zh','dG','dGt');


% Build unconditional realization
zs = fftma_perso(gen.covar, grid_gen);
ds = reshape(G * zs(:), numel(NSigma.y), numel(NSigma.x));
zhs = reshape( W * ds(:), grid_gen.ny,grid_gen.nx);

% Build condtional realization
zcs = zh + (zs - zhs);
dcs = reshape(G * zcs(:), numel(NSigma.y), numel(NSigma.x));

% Figure
figure(55);clf; c_axis=[ min(zt(:)) max(zt(:)) ];
subplot(3,3,1); imagesc(grid_gen.x, grid_gen.y, zt); caxis(c_axis); title('zt True field');
subplot(3,3,2); imagesc(grid_gen.x, grid_gen.y, d); caxis(c_axis); title('d True field');
subplot(3,3,3); imagesc(grid_gen.x, grid_gen.y, zh); caxis(c_axis); title('zh Kriging prediction based on the known averaged');
subplot(3,3,4); imagesc(grid_gen.x, grid_gen.y, zs); caxis(c_axis); title('zs Unconditional Realization');
subplot(3,3,5); imagesc(grid_gen.x, grid_gen.y, ds); caxis(c_axis); title('ds Unconditional Realization');
subplot(3,3,6); imagesc(grid_gen.x, grid_gen.y, zhs); caxis(c_axis); title('zhs Kriging prediction based on the unconditional realization');
subplot(3,3,7); imagesc(grid_gen.x, grid_gen.y, zcs); caxis(c_axis); title('zcs Conditionale Realization');
subplot(3,3,8); imagesc(grid_gen.x, grid_gen.y, dcs); caxis(c_axis); title('dcs Conditional Realization');
subplot(3,3,9); imagesc(grid_gen.x, grid_gen.y, dcs-d);  title('Error'); caxis(c_axis);

figure(66); clf
id_z = grid_gen.x < gen.covar.range0(1)*1.5;
[gamma_zt]=variogram_gridded_perso(zt);
id_d = NSigma.x < gen.covar.range0(1)*1.5;
[gamma_d]=variogram_gridded_perso(d);
subplot(2,2,1); hold on; plot(grid_gen.x(id_z)-grid_gen.x(1), 1-Cz(1,id_z)); plot(grid_gen.x(id_z)-grid_gen.x(1), gamma_zt(id_z));
subplot(2,2,2); hold on; plot(NSigma.x(id_d)-NSigma.x(1), Cd(1)- Cd(1,id_d)); plot(NSigma.x(id_d)-NSigma.x(1), gamma_d(id_d));
subplot(2,2,3); histogram(zt(:));
subplot(2,2,4); histogram(d(:));




%% Test Areal-to-point Kriging simulation
covar.model='exponential'; covar.range0=[20 20]; covar.c0=1;covar.azimuth=0;
gen.covar= kriginginitiaite(covar);
grid_gen.x=1:100;
grid_gen.y=1:100;
[grid_gen.X,grid_gen.Y] = meshgrid(grid_gen.x, grid_gen.y);
grid_gen.nxy =numel(grid_gen.X);
grid_gen.nx=numel(grid_gen.x);
grid_gen.ny=numel(grid_gen.y);
zt = fftma_perso(gen.covar, grid_gen);
NSigma.dx=5;
NSigma.dy=5;
NSigma.x = NSigma.dx/2:NSigma.dx:grid_gen.x(end);
NSigma.y = NSigma.dy/2:NSigma.dy:grid_gen.y(end);
[NSigma.X,NSigma.Y] = meshgrid(NSigma.x, NSigma.y);

% Built G
G = zeros(numel(NSigma.X),grid_gen.nxy);
for ij=1:grid_gen.nxy
     [~,id_z]=min((grid_gen.X(ij)-NSigma.X(:)).^2 + (grid_gen.Y(ij)-NSigma.Y(:)).^2);
     G(id_z,ij)=1;%/(NSigma.dx*NSigma.dy);
end
for ij = 1:numel(NSigma.X)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end

% Compute d
d = reshape(G * zt(:), numel(NSigma.y), numel(NSigma.x));

% sample zt
zt_pt = sampling_pt(grid_gen,zt,2,100);


% kriging prediction
Cz = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [grid_gen.X(:) grid_gen.Y(:)], gen.covar);
Czd = Cz * G';
Cd = G * Czd;

Czh = covardm_perso([zt_pt.x(:) zt_pt.y(:)], [zt_pt.x(:) zt_pt.y(:)], gen.covar);
Czzh = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [zt_pt.x(:) zt_pt.y(:)], gen.covar);
Czhd = Czd(zt_pt.id,:);

CCa = [ Cd Czhd' ; Czhd Czh ];
CCb = [ Czd' ; Czzh' ];

W=zeros(grid_gen.nxy,numel(NSigma.X)+zt_pt.n);
for ij=1:grid_gen.nxy
    W(ij,:) = CCa \ CCb(:,ij);
end
zh = reshape( W * [d(:) ; zt_pt.d(:)], grid_gen.ny,grid_gen.nx);



% Build unconditional realization
zs = fftma_perso(gen.covar, grid_gen);
ds = reshape(G * zs(:), numel(NSigma.y), numel(NSigma.x));
zhs = reshape( W * [ds(:) ; zs(zt_pt.id)'], grid_gen.ny,grid_gen.nx);

% Build condtional realization
zcs = zh + (zs - zhs);
dcs = reshape(G * zcs(:), numel(NSigma.y), numel(NSigma.x));


% Figure
figure(55);clf; 
c_axis=[ min(zt(:)) max(zt(:))];
subplot(3,3,1); imagesc(grid_gen.x, grid_gen.y, zt); caxis(c_axis); title('zt True field'); hold on; scatter(zt_pt.x,zt_pt.y,'xr');
subplot(3,3,2); imagesc(grid_gen.x, grid_gen.y, d); caxis(c_axis); title('d True field');
subplot(3,3,3); imagesc(grid_gen.x, grid_gen.y, zh); caxis(c_axis); title('zh Kriging prediction based on the known averaged');
subplot(3,3,4); imagesc(grid_gen.x, grid_gen.y, zs); caxis(c_axis); title('zs Unconditional Realization');
subplot(3,3,5); imagesc(grid_gen.x, grid_gen.y, ds); caxis(c_axis); title('ds Unconditional Realization');
subplot(3,3,6); imagesc(grid_gen.x, grid_gen.y, zhs); caxis(c_axis); title('zhs Kriging prediction based on the unconditional realization');
subplot(3,3,7); imagesc(grid_gen.x, grid_gen.y, zcs); caxis(c_axis); title('zcs Conditionale Realization');
subplot(3,3,8); imagesc(grid_gen.x, grid_gen.y, dcs); caxis(c_axis); title('dcs Conditional Realization');
subplot(3,3,9); imagesc(grid_gen.x, grid_gen.y, dcs-d);  title('Error'); caxis(c_axis);


figure(66); clf
id_z = grid_gen.x < gen.covar.range0(1)*1.5;
[gamma_zt]=variogram_gridded_perso(zt);
id_d = NSigma.x < gen.covar.range0(1)*1.5;
[gamma_d]=variogram_gridded_perso(d);
subplot(2,2,1); hold on; plot(grid_gen.x(id_z)-grid_gen.x(1), 1-Cz(1,id_z)); plot(grid_gen.x(id_z)-grid_gen.x(1), gamma_zt(id_z));
subplot(2,2,2); hold on; plot(NSigma.x(id_d)-NSigma.x(1), Cd(1)- Cd(1,id_d)); plot(NSigma.x(id_d)-NSigma.x(1), gamma_d(id_d));
subplot(2,2,3); histogram(zt(:));
subplot(2,2,4); histogram(d(:));

%% Reals data Areal-to-point Kriging simulation

load('result-A2PK/GEN-Run_1_2017-05-07_14-37');

% Nscore of Sigma. 
[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
NSigma.x=Sigma.x_raw; NSigma.y=Sigma.y_raw; 
[NSigma.X, NSigma.Y] = meshgrid(Sigma.x_raw, Sigma.y_raw);
d = reshape(Nscore.forward(Sigma.d_raw(:)) ,numel(NSigma.y),numel(NSigma.x));
zt = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);

% Built G
G = zeros(numel(NSigma.X),grid_gen.nxy);
for ij=1:grid_gen.nxy
     [~,id_z]=min((grid_gen.X(ij)-NSigma.X(:)).^2 + (grid_gen.Y(ij)-NSigma.Y(:)).^2);
     G(id_z,ij)=1;%/(NSigma.dx*NSigma.dy);
end
for ij = 1:numel(NSigma.X)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end


% sample zt
zt_pt = sampling_pt(grid_gen,zt,1,100);


% kriging prediction
Cz = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [grid_gen.X(:) grid_gen.Y(:)], gen.covar);
Czd = Cz * G';
Cd = G * Czd;

Czh = covardm_perso([zt_pt.x(:) zt_pt.y(:)], [zt_pt.x(:) zt_pt.y(:)], gen.covar);
Czzh = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [zt_pt.x(:) zt_pt.y(:)], gen.covar);
Czhd = Czd(zt_pt.id,:);

CCa = [ Cd Czhd' ; Czhd Czh ];
CCb = [ Czd' ; Czzh' ];

W=zeros(grid_gen.nxy,numel(NSigma.X)+zt_pt.n);
for ij=1:grid_gen.nxy
    W(ij,:) = CCa \ CCb(:,ij);
end
zh = reshape( W * [d(:) ; zt_pt.d(:)], grid_gen.ny,grid_gen.nx);



% Build unconditional realization
zs = fftma_perso(gen.covar, grid_gen);
ds = reshape(G * zs(:), numel(NSigma.y), numel(NSigma.x));
zhs = reshape( W * [ds(:) ; zs(zt_pt.id)'], grid_gen.ny,grid_gen.nx);

% Build condtional realization
zcs = zh + (zs - zhs);
dcs = reshape(G * zcs(:), numel(NSigma.y), numel(NSigma.x));


% Figure
figure(55);clf; 
c_axis=[ min(zt(:)) max(zt(:))];
subplot(3,3,1); imagesc(grid_gen.x, grid_gen.y, zt); caxis(c_axis); title('zt True field'); hold on; scatter(zt_pt.x,zt_pt.y,'xr');
subplot(3,3,2); imagesc(grid_gen.x, grid_gen.y, d); caxis(c_axis); title('d True field');
subplot(3,3,3); imagesc(grid_gen.x, grid_gen.y, zh); caxis(c_axis); title('zh Kriging prediction based on the known averaged');
subplot(3,3,4); imagesc(grid_gen.x, grid_gen.y, zs); caxis(c_axis); title('zs Unconditional Realization');
subplot(3,3,5); imagesc(grid_gen.x, grid_gen.y, ds); caxis(c_axis); title('ds Unconditional Realization');
subplot(3,3,6); imagesc(grid_gen.x, grid_gen.y, zhs); caxis(c_axis); title('zhs Kriging prediction based on the unconditional realization');
subplot(3,3,7); imagesc(grid_gen.x, grid_gen.y, zcs); caxis(c_axis); title('zcs Conditionale Realization');
subplot(3,3,8); imagesc(grid_gen.x, grid_gen.y, dcs); caxis(c_axis); title('dcs Conditional Realization');
subplot(3,3,9); imagesc(grid_gen.x, grid_gen.y, dcs-d);  title('Error'); caxis(c_axis);


figure(66); clf
id_z = grid_gen.x < gen.covar.range0(1)*1.5;
[gamma_zt]=variogram_gridded_perso(zt);
id_d = NSigma.x < gen.covar.range0(1)*1.5;
[gamma_d]=variogram_gridded_perso(d);
subplot(2,2,1); hold on; plot(grid_gen.x(id_z)-grid_gen.x(1), 1-Cz(1,id_z)); plot(grid_gen.x(id_z)-grid_gen.x(1), gamma_zt(id_z));
subplot(2,2,2); hold on; plot(NSigma.x(id_d)-NSigma.x(1), Cd(1)- Cd(1,id_d)); plot(NSigma.x(id_d)-NSigma.x(1), gamma_d(id_d));
subplot(2,2,3); histogram(zt(:));
subplot(2,2,4); histogram(d(:));



