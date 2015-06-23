function G = meas_G_grid(grid,g_true,method_G,plotit)
% MEAS_G_GRID simulate the measurement of low resolution of the grid from the g
% field (g_true). Two method are implemented (data or ERT simulation). In
% the ERT option, ERT is simulated with ERT2D package. We still need to
% implement the inverse modeling and noise to get the low resolution G
% field
% INPUT:
%       - grid:     grid of the matrix to generate (see Import_dat for more info)
%       - g_true:   variogram (see mGstat or gStat)
%       - method:   choice of method : 1. borehole | 2. random generated
%       - plotit:   1 or 0 to disply a plot or not
%
% OUTPUT:
%       - G.d:      grid measurement (G)
%       - G.std     std error associated to G
%
% Author: Raphael Nussbaumer
% date : January 2014
% need to do : add assert() for input, flexible number of input var

switch method_G
    case 1 % read data
        
        load('data.mat')
        [X,Y]=meshgrid(x,y);
        F = scatteredInterpolant(X(:),Y(:),sigma_obs(:),'linear','nearest');
        G.d = F(grid.X,grid.Y);
        F = scatteredInterpolant(X(:),Y(:),sigma_obs_err(:),'linear','nearest');
        G.err = F(grid.X,grid.Y);
        assert(~any(isnan(G.d(:))),'error');
        assert(~any(isnan(G.err(:))),'error');
        clear sigma_obs sigma_obs_err sigma_true
        
        G.std = G.err/100.*G.d;
        
        if plotit
            figure;hold on
            subplot(3,1,1); pcolor(grid.x,grid.y,g_true);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}')
            subplot(3,1,2); pcolor(grid.x,grid.y,G.d);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('G')
            subplot(3,1,3); pcolor(grid.x,grid.y,G.std);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('G_{std}')
        end
        

        
        
    case 2 % add simple noise
        noise=std(g_true(:))*.05*randn(size(g_true));
        G_noised=g_true+noise;
        
        ratio_y=5;
        ratio_x=ratio_y*12;
        
        G_upscaled = upscaling(G_noised,ratio_y,ratio_x,'arithmetique');
        G.d = downscaling(G_upscaled,ratio_y,ratio_x);
        G.std=std(G.d(:))*ones(size(noise));
        
        if plotit
            figure;hold on
            subplot(2,2,1); pcolor(grid.x,grid.y,g_true);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}')
            subplot(2,2,2); pcolor(grid.x,grid.y,G_noised);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+noise')
            subplot(2,2,3); pcolor(grid.x(1:ratio_x:end), grid.y(1:ratio_y:end), G_upscaled);shading flat;  xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+noise upscaled')
            subplot(2,2,4); pcolor(grid.x,grid.y,G.d); shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+scale upscaled and then downscaled in a grid')
        end
  
    case 3
        %error('method 3 not available yet...')
        
        script_modelling_raph(g_true)
        
        
end


G.x=grid.x;
G.y=grid.y;
G.xy=grid.xy;

end