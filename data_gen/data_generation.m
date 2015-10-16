function [grid, K_true, phi_true, sigma_true, K, sigma, Sigma, gen] = data_generation(gen)
% DATA_GENERATION is basically creating all the data required for a simulation.
% INPUT:
%       - grid{end}:         grid{end} of the matrix to generate
%           -grid{end}.nx:   nb of grid{end} cell in x
%           -grid{end}.x:    vector of the cell location [m]
%           -grid{end}.dx:   cell size [m]
%       - method:       method to get data
%           -gen:       rho_true, K_true and phi_true
%               -1:     read mat file
%               -2:     first K then phi and then g
%               -3:     first phi and then K and g
%           -samp:      K and g
%               -1:     borehole
%               -2:     random
%       - plotit:       choose to plot or not
% OUTPUT:
%       - K_true:      	Hydraulic conductivity true field, matrix (grid{end}.nx x grid{end}.ny) (data or generated)
%       - rho_true:      	Electrical conductivity true field, matrix (grid{end}.nx x grid{end}.ny)  (data or from K_true)
%       - K:          	Hydraulic conductivity at some point, structure: location (K.x, K.y) of data (K.d) (sampled from K_true)
%       - g:          	Electrical conductivity at some point, structure: location (g.x, g.y) of data (g.d) (sampled from rho_true)
%       - G:          	Electrical conductivity measured grid{end}, matrix (G.nx x G.ny) of data (G.d) (ERT inverse)
%
% Author: Raphael Nussbaumer
% date : January 2014
% need to do : add assert() for input, flexible number of input var


%% grid{end}
gen.scale.n = numel(gen.scale.x); % number of scale, ie nb of scale
grid = cell(gen.scale.n,1);

% construction of the grid at each scale
for i=1:gen.scale.n;
    grid{i}.nx=2^gen.scale.x(i)+1;
    grid{i}.ny=2^gen.scale.y(i)+1;
    grid{i}.nxy=grid{i}.nx*grid{i}.ny; % total number of cells
    
    grid{i}.dx=gen.xmax/(grid{i}.nx-1);
    grid{i}.dy=gen.ymax/(grid{i}.ny-1);
    
    grid{i}.x=linspace(0, gen.xmax, grid{i}.nx); % coordinate of cells center
    grid{i}.y=linspace(0, gen.ymax, grid{i}.ny);
    grid{i}.xy=1:grid{i}.nxy;
    
    [grid{i}.X, grid{i}.Y] = meshgrid(grid{i}.x,grid{i}.y); % matrix coordinate
end



%% 
% handle function for generating a fiel and all phsical relationship
f_new_field 	= @(grid,covar) fftma(grid.x(1),grid.dx,grid.x(end),grid.y(1),grid.dy,grid.y(end),covar);
f_Heinz         = @(phi) 10.^(6.66 *phi - 4.97); % log_10(K) = 6.66 \phi - 4.97 + noise (Heinz et al., 2003)
f_Heinz_inv     = @(K) (log10(K)+4.97)/ 6.66 ; % log_10(K) = 6.66 \phi - 4.97  + noise (Heinz et al., 2003)
f_Archie        = @(phi)43*real(phi.^1.4);  % \sigma = \sigma_W \phi ^m  + noise (Archie, 1942) where sigma_W can go up to .075, 1.2<m<1.6 
f_Archie_inv    = @(sigma) (sigma/43).^(1/1.4) ;  % \sigma = \sigma_W \phi ^m  + noise (Archie, 1942)

f_Kozeny        = @(phi,d) d^2/180*phi^3/(1-phi)^2;
f_Kozeny        = @(K,d) roots([-d^2/180/K 1 -2 1]);
f_KC            = @(phi,d10) 9810/0.001002 * phi.^3./(1-phi).^2 .* d10^2/180; % Kozeny-Carman @20°C


switch gen.method
    case 'Paolo' % Paolo's given data
        load('data.mat');
        [X,Y] = meshgrid(x,y);
        F = scatteredInterpolant(X(:),Y(:),sigma_true(:),'linear','nearest');
        rho_true= F(grid{end}.X,grid{end}.Y);
        %rho_true=interp2(x,y,sigma_true,grid{end}.X,grid{end}.Y,'nearest','extrap');
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
        K_true_n=f_new_field(grid{end},gen.covar);
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

        K = sampling_pt(grid{end},K_true,gen.samp); % 3. Simulate high-resolution point measurement of K
        g = sampling_pt(grid{end},rho_true,gen.samp); % 4. Simulate high-resolution point measurement of g
        [G, gen.G] = meas_G_grid(grid{end},rho_true,gen.G,gen.plotit); % 5. Simulate low-resolution grid{end} measurement of G

    case 'fromRho' % from phi
        field_raw   = f_new_field(grid{end},gen.covar);
        field       = ( field_raw-mean(field_raw(:)) )/ std(field_raw(:));
        phi_true    = gen.mu + gen.std*field;
        
        assert(all(phi_true(:)>0),'All phi_true are not greater than 0')
        K_true      = f_Heinz(phi_true);
        sigma_true  = f_Archie(phi_true); % archie gives conductivity, I want resisitivitiy
        rho_true    = 1000./sigma_true;
                
        K           = sampling_pt(grid{end},K_true,gen.samp); % 3. Simulate high-resolution point measurement of K
        sigma       = sampling_pt(grid{end},sigma_true,gen.samp); % 4. Simulate high-resolution point measurement of g
        rho         = sigma;
        rho.d       = 1000./sigma.d;
        % figure;imagesc(1./rho_true);axis equal;colorbar;
        [Sigma, gen] = Rho_generation(grid{end},rho_true,rho,gen); % 5. Simulate low-resolution grid{end} measurement of G
end

%% SAVING
if gen.saveit
    save(['data_gen/data/', gen.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM'), '.mat'], 'phi_true', 'K_true', 'sigma_true', 'K', 'sigma', 'Sigma','grid', 'gen')
end

fprintf('  ->  finish in %g sec\n', toc)
end