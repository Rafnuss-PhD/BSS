function S = variogram_gridded(X,grid,covar)

X=(X-mean(X(:)))./std(X(:));

% Vertical
nrbins=30;
vario_v = @(h) semivariogram1D(h,covar.c(1),covar.modele(1,3),'sph',covar.c(2));
XX=nan(grid.nx,nrbins);

for i=1:grid.nx
    temp=(X(:,i)-mean(X(:,i)))/std(X(:,i));
    Emp = variogram(grid.y',temp,'nrbins',nrbins,'plotit',false,'maxdist',15,'subsample',20000);
    XX(i,:)=Emp.val;
end
figure; subplot(1,2,1);hold on
plot(Emp.distance,mean(XX))
plot(Emp.distance,mean(XX)+std(XX))
plot(Emp.distance,mean(XX)-std(XX))
plot(Emp.distance,vario_v(Emp.distance))
legend('simulated ','theorical')
ylabel('vertical')
xlabel('m')


% Horizontal
nrbins=300;
vario_h = @(h) semivariogram1D(h,covar.c(1),covar.modele(1,2),'sph',covar.c(2));
XX=nan(grid.ny,nrbins);

for i=1:grid.ny
    temp=(X(i,:)-mean(X(i,:)))/std(X(i,:));
    Emp = variogram(grid.x',temp','nrbins',nrbins,'plotit',false,'maxdist',150,'subsample',20000);
    XX(i,:)=Emp.val;
end

subplot(1,2,2);hold on
y=mean(XX);
plot(Emp.distance(~isnan(y)),y(~isnan(y)))
plot(Emp.distance,vario_h(Emp.distance))
legend('simulated ','theorical')
ylabel('horizontal')
xlabel('m')




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
