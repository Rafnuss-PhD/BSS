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
            subplot(3,1,1); imagesc(grid.x,grid.y,g_true);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}')
            subplot(3,1,2); imagesc(grid.x,grid.y,G.d);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('G')
            subplot(3,1,3); imagesc(grid.x,grid.y,G.std);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('G_{std}')
        end
        
        
        
        
    case 2 % add simple noise

        noise=std(g_true(:))*randn(size(g_true));
        G_noised=g_true+1*noise;
        
        grid_new.nx=20;
        grid_new.ny=5;
        grid_new.x=linspace(grid.x(1),grid.x(end),grid_new.nx);
        grid_new.y=linspace(grid.y(1),grid.y(end),grid_new.ny);
        [grid_new.X, grid_new.Y]=meshgrid(grid_new.x,grid_new.y);
        dist= @(x,y) min(sqrt(sum(bsxfun(@minus,y,x).^2,2)));
        G.d=nan(size(g_true));
        for i=1:grid.nx
            for j=1:grid.ny
                [~,idx]=dist([grid.x(i) grid.y(j)],[grid_new.X(:) grid_new.Y(:)]);
                G.d(j,i)=idx;
            end
        end
        for i=1:grid_new.nx*grid_new.ny
                 G.d(G.d==i)=mean(g_true(G.d==i));
        end
        
        G.std=1.5*ones(size(G.d));
        
     
        
        if plotit
            figure;hold on
            subplot(2,2,1); imagesc(grid.x,grid.y,g_true);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}')
            subplot(2,2,2); imagesc(grid.x,grid.y,G_noised);shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+noise')
            subplot(2,2,3); imagesc(grid.x(1:ratio_x:end), grid.y(1:ratio_y:end), G_upscaled);shading flat;  xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+noise upscaled')
            subplot(2,2,4); imagesc(grid.x,grid.y,G.d); shading flat; xlabel('x[m]'); ylabel('y [m]'); title('g_{true}+scale upscaled and then downscaled in a grid')
        end
        
    case 3
        error('method 3 not available yet...')
        
        
        
end


G.x=grid.x;
G.y=grid.y;
G.xy=grid.xy;

end