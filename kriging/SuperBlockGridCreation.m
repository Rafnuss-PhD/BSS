function [k, X]=SuperBlockGridCreation(k, nx, ny, xmax, ymax, X)
%% SuperBlockGridCreation return the super block grid (k.sb)
%
% INPUT:
%
% * k       : kriging information
% * nx,ny   : number of cell (in x and y)
% * x-ymax  : regular grid max length
% * X       : Primary variable
%
% OUTPUT:
%
% * k       : kriging information with sb grid info
% * X       : Primary variable with assign sb position for all pt
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015


k.qs=[1 1; -1 1; -1 -1; 1 -1]; % used later on for orientation of the quarter windows (4 quandrant)
k.sb.nx=nx;
k.sb.ny=ny;
k.sb.x=linspace(0,xmax,k.sb.nx+1);
k.sb.dx=k.sb.x(2)-k.sb.x(1);
k.sb.x=k.sb.x(1:end-1)+k.sb.dx/2; % x will start not at 0 but at dx/2
k.sb.y=linspace(0,ymax,k.sb.ny+1);
k.sb.dy=k.sb.y(2)-k.sb.y(1);
k.sb.y=k.sb.y(1:end-1)+k.sb.dy/2;

% Creation of the superblock grid windows search 
[el_X, el_Y] = meshgrid(0:max(ceil(k.range(1)*k.wradius/k.sb.dx),ceil(k.range(2)*k.wradius/k.sb.dy)));% grid of searching windows in supergrid unit. this is a quadrant
[el_X_T, el_Y_T]=rotredtrans(el_X*k.sb.dx, el_Y*k.sb.dy, k.ang, k.range); % transforms the grid in unit
el_dist = sqrt(el_X_T.^2 + el_Y_T.^2); % distence from the point 0,0
el_X_s=el_X(el_dist<k.wradius); el_Y_s=el_Y(el_dist<k.wradius); % All point inside the windows search

% Assignment of all primary data to their super grid
X.sb_x = min([round((X.x-k.sb.x(1))/k.sb.dx +1)'; k.sb.nx*ones(1,X.n)]); % min is to avoid pb when X.x=k.sb.n
X.sb_y = min([round((X.y-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny*ones(1,X.n)]);

% Creation of a mask to select all point which belong to the windows search
% for all superblock grid cell.
k.sb.mask=false(k.sb.ny,k.sb.nx,X.n);
for i=1:k.sb.nx
    for j=1:k.sb.ny % for each cell of the supergrid...
        for u=1:length(el_X_s) %.. look at all point...
            for q=1:4 %... and assign it to the corresponding quadrant
                k.sb.mask(j,i,i+k.qs(q,1)*el_X_s(u)==X.sb_x' & j+k.qs(q,2)*el_Y_s(u)==X.sb_y')=true;
            end
        end
    end
end

% hold on; i=10; j=6;
% mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.nx+1,k.sb.ny+1))
% plot([X.x'; k.sb.x(X.sb_x)],[X.y'; k.sb.y(X.sb_y)])
% plot(X.x, X.y,'o')
% plot(X.x(k.sb.mask(j,i,:)), X.y(k.sb.mask(j,i,:)),'x','lineWidth',3)
% plot(k.sb.x(i), k.sb.y(j),'or')
% axis equal