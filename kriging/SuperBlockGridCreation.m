function [k, X]=SuperBlockGridCreation(k, xmax, ymax, X, plotit)
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
k.sb.x=linspace(0,xmax,k.sb.nx+1);
k.sb.dx=k.sb.x(2)-k.sb.x(1);
k.sb.x=k.sb.x(1:end-1)+k.sb.dx/2; % x will start not at 0 but at dx/2
k.sb.y=linspace(0,ymax,k.sb.ny+1);
k.sb.dy=k.sb.y(2)-k.sb.y(1);
k.sb.y=k.sb.y(1:end-1)+k.sb.dy/2;

% Creation of the superblock grid windows search
[el_X, el_Y] = meshgrid(0:max(ceil(k.range(1)*k.wradius/k.sb.dx),ceil(k.range(2)*k.wradius/k.sb.dy)));% grid of searching windows in supergrid unit. this is a quadrant
[el_X_T, el_Y_T]=rotredtrans(el_X*k.sb.dx, el_Y*k.sb.dy, k.rotation, k.range); % transforms the grid in unit
el_dist = sqrt(el_X_T.^2 + el_Y_T.^2); % distence from the point 0,0
k.el_X_s=el_X(el_dist<k.wradius); 
k.el_Y_s=el_Y(el_dist<k.wradius); % All point inside the windows search

% Assignment of all primary data to their super grid
X.sb_x = min([round((X.x-k.sb.x(1))/k.sb.dx +1)'; k.sb.nx*ones(1,X.n)]); % min is to avoid pb when X.x=k.sb.n
X.sb_y = min([round((X.y-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny*ones(1,X.n)]);

% Creation of a mask to select all point which belong to the windows search
% for all superblock grid cell.

k.sb.mask=false(k.sb.ny,k.sb.nx,X.n);

for i=1:k.sb.nx
    for j=1:k.sb.ny % for each cell of the supergrid...
        kt=false(X.n,1);
        for u=1:length(k.el_X_s) %.. look at all point...
            for q=1:4 %... and assign it to the corresponding quadrant
                kt(i+k.qs(q,1)*k.el_X_s(u)==X.sb_x' & j+k.qs(q,2)*k.el_Y_s(u)==X.sb_y')=true;
            end
        end
        % k.sb.mask(j,i,datasample(find(kt),min(sum(kt),nb_max),'replace',false))=true;
        k.sb.mask(j,i,kt)=true;
    end
end


if plotit
    i=5; j=6;
    windows=false(k.sb.ny,k.sb.nx);
    for u=1:length(k.el_X_s) %.. look at all point...
        for q=1:4 %... and assign it to the corresponding quadrant
            if i+k.qs(q,1)*k.el_X_s(u)<=k.sb.nx && j+k.qs(q,2)*k.el_Y_s(u)<=k.sb.ny && i+k.qs(q,1)*k.el_X_s(u)>=1 && j+k.qs(q,2)*k.el_Y_s(u)>=1% check to be inside the grid
                windows(i+k.qs(q,1)*k.el_X_s(u), j+k.qs(q,2)*k.el_Y_s(u))=true;
            end
        end
    end
    
    if plotit
        imagesc(k.sb.x,k.sb.y,windows)
        mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.nx+1,k.sb.ny+1),'EdgeColor','k','facecolor','none')
        plot(X.x, X.y,'d')
        plot(k.sb.x(i), k.sb.y(j),'or')
        plot(X.x(k.sb.mask(j,i,:)), X.y(k.sb.mask(j,i,:)),'x','lineWidth',3)
        plot([X.x'; k.sb.x(X.sb_x)],[X.y'; k.sb.y(X.sb_y)])
        axis equal tight
        legend('Super Grid','Hard data',['Center of grid ' num2str(i) ';' num2str(j)],['Hard data selected with grid ' num2str(i) ';' num2str(j)])
    end
end

