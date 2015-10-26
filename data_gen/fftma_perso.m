function field_f=fftma_perso(covar, grid)

%% Change covar if wanted
covar2 = covar;
covar2.c = 1; % [.7;0.3];
covar2.modele = gen.covar.modele(1,:);

        %% Create super grid 
grid_s.x_min = grid.x(1);
grid_s.x_max = grid.x(end)*3;
grid_s.dx    = grid.dx;
grid_s.y_min = grid.y(1);
grid_s.y_max = grid.y(end)*3;
grid_s.dy    = grid.dy;


%% Generate field

% addpath('C:\Users\rafnu\Documents\MATLAB\mGstat')
% addpath('C:\Users\rafnu\Documents\MATLAB\mGstat\misc')
% Va=[num2str(gen.covar.c(2)),' Nug(0) + ', num2str(gen.covar.c(1)),' Sph(', num2str(gen.covar.modele(1,2)), ',90,', num2str(gen.covar.modele(1,3)/gen.covar.modele(1,2)) ,')']; % V = ’sill Sph(range,rotation,anisotropy_factor)’
% field_g=fft_ma_2d(grid_s.x,grid_s.y,Va);

field_g       = fftma(grid_s.x_min,grid_s.dx,grid_s.x_max,grid_s.y_min,grid_s.dy,grid_s.y_max,covar);


%% Resample the field to initial size
field_p=field_g(grid.ny+1:2*grid.ny,grid.nx+1:2*grid.nx);


%% Adjust the field
field_f = (field_p-mean(field_p(:)))./std(field_p(:));



%% Plot
% figure;imagesc(field_f);axis equal;colorbar;
% figure; hist(field_f(:));
% variogram_gridded(field_f,grid,gen.covar)
