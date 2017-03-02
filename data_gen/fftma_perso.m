function field_f=fftma_perso(covar, grid)

%% Create super grid 
grid_s.x_min = grid.x(1);
grid_s.x_max = grid.x(end)*3;
grid_s.y_min = grid.y(1);
grid_s.y_max = grid.y(end)*3;

if ~isfield(grid, 'dx')
    grid_s.dx    = grid.x(2)-grid.x(1);
    grid_s.dy    = grid.y(2)-grid.y(1);
    
else
    grid_s.dx    = grid.dx;
    grid_s.dy    = grid.dy;
end

if ~isfield(grid, 'nx')
    grid.nx = numel(grid.x);
    grid.ny = numel(grid.y);
end



%% Generate field

% addpath('C:\Users\rafnu\Documents\MATLAB\mGstat')
% addpath('C:\Users\rafnu\Documents\MATLAB\mGstat\misc')
% Va=[num2str(gen.covar.c(2)),' Nug(0) + ', num2str(gen.covar.c(1)),' Sph(', num2str(gen.covar.modele(1,2)), ',90,', num2str(gen.covar.modele(1,3)/gen.covar.modele(1,2)) ,')']; % V = �sill Sph(range,rotation,anisotropy_factor)�
% field_g=fft_ma_2d(grid_s.x,grid_s.y,Va);
field_g       = fftma(grid_s.x_min,grid_s.dx,grid_s.x_max,grid_s.y_min,grid_s.dy,grid_s.y_max,covar);


%% Resample the field to initial size
field_p=field_g(grid.ny+1:2*grid.ny,grid.nx+1:2*grid.nx);


%% Adjust the field
field_f = (field_p-mean(field_p(:)))./std(field_p(:));



%% Plot

% figure;imagesc(field_f);axis equal;colorbar;
% figure; hist(field_f(:));
% 
% myfun = @(x,h) semivariogram1D(h,1,x,'sph',0);
% 
% [gamma_x, gamma_y] = variogram_gridded_perso(field_f);
% figure; subplot(1,2,1); hold on;
% id= grid.x<covar.modele(1,2);
% plot(grid.x(id),gamma_x(id),'linewidth',2)
% plot(grid.x(id),myfun(covar.modele(1,2),grid.x(id)),'linewidth',2)
% subplot(1,2,2);hold on
% id= grid.y<covar.modele(1,3);
% plot(grid.y(id),gamma_y(id),'linewidth',2)
% plot(grid.y(id),myfun(covar.modele(1,3),grid.y(id)),'linewidth',2)