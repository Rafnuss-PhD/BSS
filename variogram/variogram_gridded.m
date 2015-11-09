function [x,y] = variogram_gridded(X,grid,range,nrbins,subsample_var,subsample_grid)

X=(X-mean(X(:)))./std(X(:));



% Horizontal
xx = nan(min([subsample_grid(1),grid.ny]),nrbins(1));
for i=ceil(linspace(1,grid.ny, min([subsample_grid(1),grid.ny]) ))
    temp=(X(i,:)-mean(X(i,:)))/std(X(i,:));
    Emp = variogram(grid.x',temp','nrbins',nrbins(1),'plotit',false,'maxdist',range(1)*1.3,'subsample',subsample_var(1));
    xx(i,:)=Emp.val;
end
x.val = nanmean(xx);
x.dist = Emp.distance;

% Vertical
yy=nan(min([subsample_grid(2),grid.nx]),nrbins(2));

for i=1:ceil(linspace(1,grid.nx, min([subsample_grid(2),grid.nx]) ))
    temp=(X(:,i)-mean(X(:,i)))/std(X(:,i));
    Emp = variogram(grid.y',temp,'nrbins',nrbins(2),'plotit',false,'maxdist',range(2)*1.3,'subsample',subsample_var(2));
    yy(i,:)=Emp.val;
end
y.val = nanmean(yy);
y.dist = Emp.distance;


% vario_h = @(h) semivariogram1D(h,c,range(1),'sph',covar.c(2));
% vario_v = @(h) semivariogram1D(h,c,range(2),'sph',covar.c(2));

% [ny,nx]=size(Z);
% x = (1:nx)*dx;
% y = (1:ny)*dy;
% 
% % Horizontal
% gamma_x=nan(nx,ny*(nx-1));
% for d=1:nx
%     u=1;
%     for j=1:ny
%         for i=1:(nx-d)
%             gamma_x(d,u) = ( Z(j,i)-Z(j,i+d) )^2;
%             u=u+1;
%         end
%     end
% end
% 
% % Vertical
% gamma_y=nan(ny,nx*(ny-1));
% for d=1:ny
%     u=1;
%     for i=1:nx
%         for j=1:(ny-d)
%             gamma_y(d,u) = ( Z(j,i)-Z(j+d,i) )^2;
%             u=u+1;
%         end
%     end
% end
% 
% 
% 
% figure; hold on
% %h=repmat(x(1:end-1)',1,ny*(nx-1));
% %plot(h(:),gamma(:),'x')
% subplot(2,1,1)
% plot(x',nanmean(gamma_x,2),'x')
% subplot(2,1,2)
% plot(y',nanmean(gamma_y,2),'x')
end
