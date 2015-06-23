function flow(grid,K,plotit)
%% FLOW : velocity and tracers
% FLOW is computing steady state flow in the K-field and function of
% repartion
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch) 
% * *Date:* 29.01.2015

addpath('./flow/')

% Setting basic parameters
P0=0; P1=1;
viscosite=1; phi=.5;
Nparticles=100;

%% STEADY
% Compute the steady-state velocity field
[~,ux,uy] = steady(K,P0,P1,viscosite,phi,grid.dx,grid.dy);




%% TRACER
% Compute streamline and tracer experiment
if P0 > P1
    xo = grid.dx;
else
    xo = (grid.nx-2)*grid.dx;
end
time=zeros(1,Nparticles);ti=tic;
for u = 1:Nparticles
    stop = 0;
    x = [];
    while isempty(x) && stop < 20
        yo = min(grid.y)+range(grid.y)*rand(1);
        [x,~,tof] = addstreamlineFD_Flowv2(xo,yo,grid.nx,grid.ny,grid.dx,grid.dy,ux,uy,P0,P1);
        stop = stop+1;
    end
    clear x y
    disp(['indice ',num2str(u),' over ' ,num2str(Nparticles), ' (', num2str(u/Nparticles*100),'%) ','  time  ',num2str(toc(ti)) ,' s (expected: ',num2str(toc(ti)*Nparticles/u/60), ' min - and remaining: ',num2str(toc(ti)*Nparticles/u/60-toc(ti)/60),')'])
    time(u) = tof;
end
[temps,nombre] = tracertimes(time);


%% GRID
% Interpolate the velocity from the cell interface to the center of each
% cells.

ux_x=grid.x(1:end-1)+grid.dx/2;
ux_y=grid.y;
uy_x=grid.x;
uy_y=grid.y(1:end-1)+grid.dy/2;

% Generation of the 2D grid for both ux and uy
[ux_X,ux_Y] = meshgrid(ux_x,ux_y);
[uy_X,uy_Y] = meshgrid(uy_x,uy_y);

% Uinterpolate at the usual grid ux and uy

[X,Y] = meshgrid(grid.x,grid.y);
UX = interp2(ux_X,ux_Y,ux,X,Y);
UY = interp2(uy_X,uy_Y,uy,X,Y);

% Define start point for streamline
starty = grid.y;
startx = grid.x(end-2)*ones(grid.ny,1);



%% PLOT
if plotit
    subplot(2,1,1); hold on;
    pcolor(grid.x,grid.y,log10(K)); shading flat;
    streamline(X,Y,UX,UY,startx,starty)
    xlabel('x[m]'); ylabel('y [m]'); title('Hydraulic Conductifity and streamline'); colorbar;
    
    subplot(2,1,2)
    plot(temps,nombre)
end