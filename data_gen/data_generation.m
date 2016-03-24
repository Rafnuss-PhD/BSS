function [filename, grid_gen, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen)
% DATA_GENERATION is basically creating all the data required for a simulation.
% INPUT:
% GRID
%   * xmax:         length of the grid_gen of x [unit]
%   * ymax:         length of the grid_gen of y [unit]
%   * sx:          scale level. number of cell is 2^gen.scale.x(i)+1
%   * sy:          scale level. number of cell is 2^gen.scale.y(i)+1
% STRUCTURE
%   * method:       method of generation: 'fromRho', ''
%   * samp:         Method of sampling of K and g | 1: borehole, 2:random. For fromK or from Rho only
%   * samp_n:       Number of well or number of point
%   * covar:
%       * modele    covariance structure
%       * c          
%   * mu            parameter of the first field. 
%   * std           
% RHO R2
%   * Rho.grid_gen.nx           = 300;
%   * Rho.grid_gen.ny           = 60; % log-spaced grid_gen.
%   * Rho.elec.spacing      = 2; % in grid_gen spacing unit.
%   * Rho.elec.config_max   = 6000; % number of configuration of electrode maximal 
%   * Rho.dmin.res_matrix   = 1; % resolution matrix: 1-'sensitivity' matrix, 2-true resolution matrix or 0-none
%   * Rho.dmin.tolerance    = 1;
%   * Rho.method          = 'R2';
% OTHER
%   * plotit              = false;      % display graphic or not (you can still display later with |script_plot.m|)
%   * saveit              = true;       % save the generated file or not, this will be turn off if mehod Paolo or filename are selected
%   * name                = 'Small_range';
%   * seed                = 123456;
% OUTPUT:
%       - K_true:      	Hydraulic conductivity true field, matrix (grid_gen.nx x grid_gen.ny) (data or generated)
%       - rho_true:      	Electrical conductivity true field, matrix (grid_gen.nx x grid_gen.ny)  (data or from K_true)
%       - K:          	Hydraulic conductivity at some point, structure: location (K.x, K.y) of data (K.d) (sampled from K_true)
%       - g:          	Electrical conductivity at some point, structure: location (g.x, g.y) of data (g.d) (sampled from rho_true)
%       - G:          	Electrical conductivity measured grid_gen, matrix (G.nx x G.ny) of data (G.d) (ERT inverse)
%
% Author: Raphael Nussbaumer

%% * *INPUT CEHCKING*
assert(isfield(gen, 'xmax'))
assert(isfield(gen, 'ymax'))
assert(isfield(gen, 'sx'))
assert(isfield(gen, 'sy'))
if ~isfield(gen, 'method'); gen.method = 'Random'; end
if ~isfield(gen, 'samp'); gen.samp = 2; end
if ~isfield(gen, 'samp_n'); gen.samp_n = 1/100 * (2^gen.sx+1)*(2^gen.sy+1); end
if ~isfield(gen, 'mu'); gen.mu = 0; end
if ~isfield(gen, 'std'); gen.std = 1; end
% Inversion for secondary variable
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'method'); gen.Rho.method = 'R2'; end
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'grid_gen') || ~isfield(gen.Rho.grid_gen, 'nx'); gen.Rho.grid_gen.nx = 300; end
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'grid_gen') || ~isfield(gen.Rho.grid_gen, 'ny'); gen.Rho.grid_gen.nx = 60; end
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'elec') || ~isfield(gen.Rho.elec, 'spacing'); gen.Rho.elec.spacing = 2; end
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'elec') || ~isfield(gen.Rho.elec, 'config_max'); gen.Rho.elec.config_max = 6000; end
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'dmin') || ~isfield(gen.Rho.dmin, 'res_matrix'); gen.Rho.dmin.res_matrix = 1; end
if ~isfield(gen, 'Rho') || ~isfield(gen.Rho, 'dmin') || ~isfield(gen.Rho.dmin, 'tolerance'); gen.Rho.dmin.tolerance = 1; end
% Other
if ~isfield(gen, 'plotit'); gen.plotit = 0; end
if ~isfield(gen, 'saveit'); gen.saveit = 1; end
if ~isfield(gen, 'name'); gen.name = ''; end
if ~isfield(gen, 'seed'); gen.seed = rand; end

tic

%% * 2. *construction of the grid_gen*
grid_gen.sx=gen.sx;
grid_gen.sy=gen.sy;
grid_gen.nx=2^gen.sx+1;
grid_gen.ny=2^gen.sy+1;
grid_gen.nxy=grid_gen.nx*grid_gen.ny; % total number of cells

grid_gen.dx=gen.xmax/(grid_gen.nx-1);
grid_gen.dy=gen.ymax/(grid_gen.ny-1);

grid_gen.x=linspace(0, gen.xmax, grid_gen.nx); % coordinate of cells center
grid_gen.y=linspace(0, gen.ymax, grid_gen.ny);
grid_gen.xy=1:grid_gen.nxy;

[grid_gen.X, grid_gen.Y] = meshgrid(grid_gen.x,grid_gen.y); % matrix coordinate

%% * 2. *handle function for generating a fiel and all phsical relationship*

f_Heinz         = @(phi) 10.^(6.66 *phi - 4.97); % log_10(K) = 6.66 \phi - 4.97 + noise (Heinz et al., 2003)
f_Heinz_inv     = @(K) (log10(K)+4.97)/ 6.66 ; % log_10(K) = 6.66 \phi - 4.97  + noise (Heinz et al., 2003)
f_Archie        = @(phi)43*real(phi.^1.4);  % \sigma = \sigma_W \phi ^m  + noise (Archie, 1942) where sigma_W can go up to .075, 1.2<m<1.6 
f_Archie_inv    = @(sigma) (sigma/43).^(1/1.4) ;  % \sigma = \sigma_W \phi ^m  + noise (Archie, 1942)

f_Kozeny        = @(phi,d) d^2/180*phi^3/(1-phi)^2;
f_Kozeny        = @(K,d) roots([-d^2/180/K 1 -2 1]);
f_KC            = @(phi,d10) 9810/0.001002 * phi.^3./(1-phi).^2 .* d10^2/180; % Kozeny-Carman @20°C

%% * 3. *Generate field*
switch gen.method
    case 'Paolo' % Paolo's given data
        load('data.mat');
        [X,Y] = meshgrid_gen(x,y);
        F = scatteredInterpolant(X(:),Y(:),sigma_true(:),'linear','nearest');
        rho_true= F(grid_gen.X,grid_gen.Y);
        %rho_true=interp2(x,y,sigma_true,grid_gen.X,grid_gen.Y,'nearest','extrap');
        assert(~any(isnan(rho_true(:))),'error')
        % figure;hold on;ksdensity(rho_true(:));ksdensity(sigma_true(:))
        % figure;subplot(2,1,1);imagesc(rho_true);subplot(2,1,2);imagesc(sigma_true);
        clear sigma_obs sigma_obs_err sigma_true
        
        phi_true = f_Archie_inv(rho_true);
        K_true =  f_Heinz(phi_true);
        % K_true=nan(size(rho_true)); warning('K_true not generated. read from file unavailable')
        % phi_true=nan(size(rho_true)); warning('K_true not generated. read from file unavailable')
        info.G=1;
        gen.saveit =false; % turn save off by default
        gen.G = 'Paolo';
               
    case 'fromK' % from K
        K_true_n=f_new_field(grid_gen,gen.covar);
        % figure;imagesc(K_true_n)
        K_true= 10.^(gen.mu + K_true_n * sqrt(gen.std)); % put it as a log normal dist.
        assert(all(K_true(:)>0),'All K_true are not greater than 0')
        if any(K_true(:)<10^-4.96)
            warning(['Heinz does not work for ' num2str(sum(K_true(:)<10^-4.96)) ' data'])
            K_true(K_true<10^-4.96)=10^-4.96;
        end
        phi_true=f_Heinz_inv(K_true);
        assert(all(phi_true(:)>0),'All phi_true are not greater than 0')
        s_true = f_Archie(phi_true);
        rho_true          = 1./s_true;

        K = sampling_pt(grid_gen,K_true,gen.samp); % 3. Simulate high-resolution point measurement of K
        g = sampling_pt(grid_gen,rho_true,gen.samp); % 4. Simulate high-resolution point measurement of g
        [G, gen.G] = meas_G_grid_gen(grid_gen,rho_true,gen.G,gen.plotit); % 5. Simulate low-resolution grid_gen measurement of G
        
    case 'Normal-Random'
        field = fftma_perso(gen.covar, grid_gen);

        sigma_true      =  gen.mu + gen.std*field;
        
        sigma           = sampling_pt(grid_gen,sigma_true,gen.samp,gen.samp_n);
        phi_true        = NaN;
        K_true          = NaN;
        K               = NaN;
        Sigma.d         = sigma_true;
        Sigma.x         = grid_gen.x;
        Sigma.y         = grid_gen.y;
        Sigma.std       = 1+zeros(size(Sigma.d));
        
    case 'Log-Random'
        sigma_true      =  10.^(gen.mu + gen.std*fftma_perso(gen.covar, grid_gen));
        % figure;imagesc(sigma_true);axis equal;colorbar;
        sigma           = sampling_pt(grid_gen,sigma_true,gen.samp,gen.samp_n);
        phi_true        = NaN;
        K_true          = NaN;
        K               = NaN;
        Sigma.d         = sigma_true;
        Sigma.x         = grid_gen.x;
        Sigma.y         = grid_gen.y;
        Sigma.std       = 1+zeros(size(Sigma.d));
        
    case 'fromRho' % from phi
        phi_true    = gen.mu + gen.std*fftma_perso(gen.covar, grid_gen);
        
        assert(all(phi_true(:)>0),'All phi_true are not greater than 0')
        K_true      = f_Heinz(phi_true);
        sigma_true  = f_Archie(phi_true); % archie gives conductivity, I want resisitivitiy
        rho_true    = 1000./sigma_true;
                
        K           = sampling_pt(grid_gen,K_true,gen.samp,gen.samp_n); % 3. Simulate high-resolution point measurement of K
        sigma       = sampling_pt(grid_gen,sigma_true,gen.samp,gen.samp_n); % 4. Simulate high-resolution point measurement of g
        rho         = sigma;
        rho.d       = 1000./sigma.d;
        
        figure(1);clf; subplot(2,1,1);imagesc(grid_gen.x,grid_gen.y,sigma_true);colorbar; hold on; scatter(sigma.x,sigma.y,sigma.d)
        subplot(2,1,2); hold on;[f,x]=hist(sigma_true(:)); plot(x,f/trapz(x,f)); [f,x]=hist(sigma.d(:)); plot(x,f/trapz(x,f));
        
        [gamma_x, gamma_y] = variogram_gridded_perso(sigma_true); myfun = @(x,h) semivariogram1D(h,1,x,gen.covar.modele(1),0);
        figure(2); clf; subplot(2,1,1); hold on; plot(grid_gen.x(1:end/2),gamma_x(1:end/2)); plot([gen.covar.modele(1,2) gen.covar.modele(1,2)],[0 1])
        plot(grid_gen.x(1:end/2),myfun(gen.covar.modele(1,2),grid_gen.x(1:end/2)),'linewidth',2)
        subplot(2,1,2); hold on;plot(grid_gen.y(1:end/2),gamma_y(1:end/2)); plot([gen.covar.modele(1,3) gen.covar.modele(1,3)],[0 1])
        plot(grid_gen.y(1:end/2),myfun(gen.covar.modele(1,3),grid_gen.y(1:end/2)),'linewidth',2)
        
        
        keyboard;close all;
        
        [Sigma, gen] = Rho_generation(grid_gen,rho_true,gen); % 5. Simulate low-resolution grid_gen measurement of G
        
        Sigma        = Sens2std(Sigma,grid_gen,sigma);
    otherwise
        error('method not define.')
end

%% * 4.*SAVING*

if gen.saveit
    filename = ['data_gen/data/GEN-', gen.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM'), '.mat'];
    save(filename, 'phi_true', 'K_true', 'sigma_true', 'K', 'sigma', 'Sigma','grid_gen', 'gen')
end

fprintf('  ->  finish in %g sec\n', toc)
end