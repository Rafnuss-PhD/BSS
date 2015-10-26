
%% R2  
% Electrode configuration
figure;
plot(gen.Rho.elec.data,1:gen.Rho.elec.confirho_n,'x')
set(gca, 'YDir', 'reverse'); set(gca,'xtick',0:gen.Rho.elec.n)
xlabel('electrode position'); ylabel('Configuration')
legend('A','B','M','N');


% Pseudo-section: electrode position
figure; hold on
plot(1:gen.Rho.elec.n,zeros(1,gen.Rho.elec.n),'x'); axis equal
scatter(gen.Rho.elec.pos_ini(:,1),gen.Rho.elec.pos_ini(:,2),[],gen.Rho.elec.k_ini); 
plot(gen.Rho.elec.pos(:,1),gen.Rho.elec.pos(:,2),'or')
xlabel('electrode position');ylabel('depth')
set(gca, 'YDir', 'reverse');


% Apparent Resisitivity
figure;
subplot(3,1,1);
hold on; title('Observed Apparent Resistivity')
contourf(gen.Rho.grid.x,gen.Rho.grid.y, gen.Rho.f.output.pseudo_interp,'EdgeColor', 'none');
scatter(gen.Rho.elec.pseudo_x, gen.Rho.elec.pseudo_y,[],gen.Rho.f.output.pseudo,'filled','MarkerEdgeColor','k')
set(gca,'Ydir','reverse'); axis equal; axis tight;colorbar;
subplot(3,1,2);
hold on; title('Calculated Apparent Resistivity')
contourf(gen.Rho.grid.x,gen.Rho.grid.y, gen.Rho.i.output.pseudo_interp,'EdgeColor', 'none')
scatter(gen.Rho.elec.pseudo_x, gen.Rho.elec.pseudo_y,[],gen.Rho.i.output.pseudo,'filled','MarkerEdgeColor','k')
set(gca,'Ydir','reverse'); axis equal; axis tight;colorbar;
subplot(3,1,3);
hold on; title('Error in Apparent Resistivity')
contourf(gen.Rho.grid.x,gen.Rho.grid.y, gen.Rho.i.output.err_interp,'EdgeColor', 'none')
scatter(gen.Rho.elec.pseudo_x, gen.Rho.elec.pseudo_y,[],gen.Rho.i.output.err,'filled','MarkerEdgeColor','k')
set(gca,'Ydir','reverse'); axis equal; axis tight;colorbar;


%
figure;
caxis_lim = [min(gen.Rho.f.output.res(:)) max(gen.Rho.f.output.res(:))];
subplot(4,1,1);
surf(gen.Rho.grid.x,gen.Rho.grid.y,gen.Rho.f.output.res,'EdgeColor', 'none');view(2);
set(gca,'Ydir','reverse'); caxis(caxis_lim); axis tight;colorbar;% axis equal;
xlabel('x[m]'); ylabel('y [m]'); title('Upscale resisitivity rho_{true,upscaled}');
subplot(4,1,2);
surf(gen.Rho.grid.x,gen.Rho.grid.y,gen.Rho.i.output.res,'EdgeColor', 'none');view(2);
set(gca,'Ydir','reverse'); caxis(caxis_lim); axis tight; colorbar;% axis equal;
xlabel('x[m]'); ylabel('y [m]'); title('Inversed resisitivity G');
subplot(4,1,3);
surf(gen.Rho.grid.x,gen.Rho.grid.y,(gen.Rho.i.output.res-gen.Rho.f.output.res),'EdgeColor', 'none');view(2);
set(gca,'Ydir','reverse'); axis tight; colorbar;% axis equal;
xlabel('x[m]'); ylabel('y [m]'); title('Error in resistivity resisitivity rho_{true,upscaled}-G');
subplot(4,1,4);
surf(grid{end}.x, grid{end}.y, Rho.std,'EdgeColor', 'none');view(2);
set(gca,'Ydir','reverse');  axis tight; colorbar; % axis equal;
xlabel('x[m]'); ylabel('y [m]'); title('INTERPOLATED Standard deviation of resistivity Rho.std');


% Well data
figure;

title('Observed and inverted electrical resisitivity at well')
[uniq_x,ia,ic] = unique(sigma.x);  %C = A(ia) and A = C(ic).
n_well = numel(uniq_x);
Rho_d_well = Rho.d(ismember(sigma.y,grid{end}.Y)&ismember(sigma.x,grid{end}.X));
Rho_sen_well = Rho.sen(ismember(sigma.y,grid{end}.Y)&ismember(sigma.x,grid{end}.X));
for i_well=1:n_well 
    subplot(1,n_well,i_well); hold on;
    plot(Rho_d_well(ic==i_well), sigma.y(ic==i_well))
    plot(sigma.d(ic==i_well), sigma.y(ic==i_well))
    xlim([floor(min(sigma.d)) ceil(max(sigma.d))])
    ylabel('y [m]'); xlabel(['Electrical Conductivity at ' num2str(uniq_x(i_well)) 'm']); set(gca,'Ydir','reverse');
end

rho_err=Rho_d_well-sigma.d;
figure;
plot(Rho_sen_well,rho_err,'x');
xlabel('Normalized sensitivity');ylabel('Error between well measurement and inverted image')


%% Generated Data  

figure;
caxis_lim = [min([Sigma.d(:); sigma_true(:)]) max([Sigma.d(:); sigma_true(:)]) ];
subplot(4,1,1); 
imagesc(grid{end}.x, grid{end}.y, sigma_true); 
title('True Porosity \rho_{true}'); 
xlabel('x[m]'); ylabel('y [m]'); colorbar;set(gca,'Ydir','reverse');
subplot(4,1,2); hold on; 
imagesc(grid{end}.x, grid{end}.y, log10(K_true)); plot(K.x, K.y, 'or'); 
title('Log True Hydraulic Conudctivity K_{true} and sampled point location K'); 
xlabel('x[m]'); ylabel('y [m]'); colorbar; set(gca,'Ydir','reverse'); axis tight
subplot(4,1,3);
imagesc(grid{end}.x, grid{end}.y, sigma_true); hold on; plot(sigma.x, sigma.y, 'or'); 
title('True Electrical Conductivity \sigma_{true} and sampled point location g');
xlabel('x[m]'); ylabel('y [m]');  colorbar; set(gca,'Ydir','reverse'); caxis(caxis_lim);
subplot(4,1,4); 
imagesc(grid{end}.x, grid{end}.y, Sigma.d); 
title('Inverted Electrical Conductivity \Sigma_{true}'); 
xlabel('x[m]'); ylabel('y [m]');colorbar; set(gca,'Ydir','reverse'); caxis(caxis_lim);


figure;
subplot(4,1,1); 
pcolor(grid{end}.x,grid{end}.y,Sigma.std); shading flat; 
xlabel('x[m]'); ylabel('y [m]'); colorbar; set(gca,'Ydir','reverse');
title('Electrical Conductivity Tomography error Rho_{std}'); 
subplot(4,1,2); 
imagesc(grid{end}.x,grid{end}.y,Sigma.d); 
xlabel('x[m]'); ylabel('y [m]');  colorbar;set(gca,'Ydir','reverse'); caxis(caxis_lim);
title('Electrical Conductivity Tomography G');
subplot(4,1,3); 
imagesc(grid{end}.x,grid{end}.y,sigma_true);
xlabel('x[m]'); ylabel('y [m]'); colorbar;set(gca,'Ydir','reverse'); caxis(caxis_lim);
title('True Electrical Conductivity rho_{true}');
subplot(4,1,4); 
imagesc(grid{end}.x,grid{end}.y,sigma_true-Sigma.d); 
xlabel('x[m]'); ylabel('y [m]'); title('rho_{true}-G');colorbar;set(gca,'Ydir','reverse');

% Histogram
figure;
nbins_all = grid{end}.nxy/1000;
nbins_sampled = sigma.n/20;
subplot(2,2,1);
[f,x]=hist(sigma_true(:),nbins_all); plot(x,f/trapz(x,f));
title(sprintf('True porosity \\rho_{true} | \\mu=%.2f, \\sigma=%.2f',mean(sigma_true(:)), std(sigma_true(:))))
subplot(2,2,2); hold on
[f,x]=hist(log10(K_true(:)),nbins_all); plot(x,f/trapz(x,f));
[f,x]=hist(log10(K.d(:)),nbins_sampled); plot(x,f/trapz(x,f));
legend('K_{true}','sampled K')
title(sprintf('LOG permeability| \\mu=%.2f, \\sigma=%.2f',mean(log10(K_true(:))), std(log10(K_true(:)))))
subplot(2,2,3); hold on
[f,x]=hist(sigma_true(:),nbins_all); plot(x,f/trapz(x,f));
[f,x]=hist(sigma.d(:),nbins_sampled); plot(x,f/trapz(x,f));
legend('rho_{true}','sampled g')
title(sprintf('Resisitivity | \\mu=%.2f, \\sigma=%.2f',mean(sigma_true(:)), std(sigma_true(:))))
subplot(2,2,4); 
[f,x]=hist(Sigma.d(:),nbins_all); plot(x,f/trapz(x,f));
title(sprintf('Inversed resisitivity Rho_{true} | \\mu=%.2f, \\sigma=%.2f',mean(Sigma.d(:)), std(Sigma.d(:))))


%%  OUPTUT  

% 
% %% Kriging
% figure(1); clf
% imagesc(Y.x,Y.y,Y.m_ns{1}); axis tight; lim=axis;hold on
% 
% t=-pi:0.01:pi;
% x=Y.x(Y.pt.x)+k.range(1)*cos(t);
% y=Y.y(Y.pt.y)+k.range(2)*sin(t);
% plot(x,y,'--r')
% axis(lim)
% % scatter(Y.X(~isnan(Y.m_ns{1})),Y.Y(~isnan(Y.m_ns{1})),[],Y.m_ns{1}(~isnan(Y.m_ns{1})),'d','filled','MarkerEdgeColor','k'); axis tight
% % 
% % scatter(X.x, X.y,'o','k') % well sampling
% 
% lambda_c= 36+50.*(k0.lambda-min(k0.lambda))./range(k0.lambda);
% XY_ns = [X.d_ns; Y.m_ns{1}(~isnan(Y.m_ns{1}))];
% 
% scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns(k0.mask),'o','filled','MarkerEdgeColor','k'); axis tight
% 
% 
% scatter(Y.x(Y.pt.x),Y.y(Y.pt.y),100,k0.lambda'* XY_ns(k0.mask),'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
% 
% 
% % scatter(X.x(k0.mask), X.y(k0.mask),[],lambda_c(1:sum(k0.sb_mask)),'o','filled')
% % scatter(Y.X(k0.mask), Y.Y(k0.mask),[],lambda_c(sum(k0.sb_mask)+1:end),'d','filled')
% % 
% % scatter(X.x(k0.sb_mask), X.y(k0.sb_mask),[],lambda_c(1:sum(k0.sb_mask)),'o','filled')
% % scatter(Y.X(k0.ss_mask), Y.Y(k0.ss_mask),[],lambda_c(sum(k0.sb_mask)+1:end),'d','filled')
% 
% 
% legend('pt estimated','hard data selected','soft data selected')
% axis tight;


%% BSGS
% Multi-grid
sub_n = min([parm.n_realisation,5])+1;
if sub_n<4, kk=sub_n;mm=1;
else kk=ceil(sub_n/2); mm=2; 
end

figure;
caxis_limm = [min([X_true(:);Y{end}.m{end}(:)]) max([X_true(:);Y{end}.m{end}(:)])];
subplot(kk,mm,1); hold on; 
pcolor(grid{end}.x,grid{end}.y,X_true);
plot(X.x,X.y,'or')
shading flat;colorbar; caxis(caxis_limm); axis tight;set(gca,'Ydir','reverse'); xlabel('x[m]'); ylabel('y[m]');
title([ 'True 1/',parm.unit ' X_{true}: \mu=' num2str(mean2(X_true)) ' | \sigma='   num2str(std2(X_true)) ' and sampled g: \mu=' num2str(mean(X.d)) ' | \sigma='   num2str(std(X.d))]);
for i_sim=2:mm*kk
    subplot(kk,mm,i_sim); 
    pcolor(grid{parm.scale(end)}.x,grid{parm.scale(end)}.y,Y{end}.m{i_sim-1});
    shading flat;colorbar;caxis(caxis_limm);axis tight;set(gca,'Ydir','reverse'); xlabel('x[m]');ylabel('y[m]');
    title(['Simmulated 1/', parm.unit, ' gG: \mu=' num2str(mean2(Y{end}.m{i_sim-1})) ' | \sigma=' num2str(std2(Y{end}.m{i_sim-1})) ]);
end


% Histogramm
figure; 
nbins=50;
subplot(2,1,1);hold on; title('Result')
[f,x]=hist(X_true(:)); plot(x,f/trapz(x,f),'linewidth',2);
[f,x]=hist(X.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
[f,x]=hist(Z.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
for i_sim=1:mm*kk-1
    [f,x]=hist(Y{end}.m{i_sim}(:),nbins); plot(x,f/trapz(x,f));
end
legend('True X', 'Sampled X', 'Z', 'realisation')

subplot(2,1,2);hold on;title('Normal Space')
[f,x]=hist(Nscore.forward(X_true(:))); plot(x,f/trapz(x,f),'linewidth',2);
[f,x]=hist(Nscore.forward(X.d(:))); plot(x,f/trapz(x,f),'linewidth',2);
[f,x]=hist(Nscore.forward(Z.d(:))); plot(x,f/trapz(x,f),'linewidth',2);
for i_sim=1:mm*kk-1
    [f,x]=hist(Nscore.forward(Y{end}.m{i_sim}(:))); plot(x,f/trapz(x,f));
end
legend('True X', 'Sampled X', 'Z', 'realisation')


%% Variogram  
nrbins=30;
myfun = @(x,h) semivariogram1D(h,1,x,'sph',0);
val_m=nan(grid{parm.scale(end)}.nx,nrbins);
val_true=nan(grid{parm.scale(end)}.nx,nrbins);

for i=1:grid{parm.scale(end)}.nx
    temp=(sigma_true(:,i)-mean(sigma_true(:,i)))/std(sigma_true(:,i));
    Emp = variogram(grid{end}.y',temp,'nrbins',nrbins,'plotit',false,'maxdist',15,'subsample',20000);
    val_true(i,:)=Emp.val;

    temp=(Y{end}.m{end}(:,i)-mean(Y{end}.m{end}(:,i)))/std(Y{end}.m{end}(:,i));
    Emp = variogram(grid{parm.scale(end)}.y',temp,'nrbins',nrbins,'plotit',false,'maxdist',15,'subsample',20000);
    val_m(i,:)=Emp.val;
end
figure; subplot(1,2,1);hold on
plot(Emp.distance,mean(val_true))
plot(Emp.distance,mean(val_m))
plot(Emp.distance,myfun(5,Emp.distance))
legend('true conductivity','simulated conductivity','theorical equation')
ylabel('vertical')
xlabel('m')

nrbins=300;
val_m=nan(grid{parm.scale(end)}.ny,nrbins);
val_true=nan(grid{parm.scale(end)}.ny,nrbins);

for i=1:grid{parm.scale(end)}.ny
    temp=(sigma_true(i,:)-mean(sigma_true(i,:)))/std(sigma_true(i,:));
    Emp = variogram(grid{end}.x',temp','nrbins',nrbins,'plotit',false,'maxdist',150,'subsample',20000);
    val_true(i,:)=Emp.val;
    
    temp=(Y{end}.m{end}(i,:)-mean(Y{end}.m{end}(i,:)))/std(Y{end}.m{end}(i,:));
    Emp = variogram(grid{parm.scale(end)}.x',temp','nrbins',nrbins,'plotit',false,'maxdist',150,'subsample',20000);
    val_m(i,:)=Emp.val;
end

subplot(1,2,2);hold on
y=mean(val_true);
plot(Emp.distance(~isnan(y)),y(~isnan(y)))
y=mean(val_m);
plot(Emp.distance(~isnan(y)),y(~isnan(y)))
plot(Emp.distance,myfun(75,Emp.distance))
legend('true conductivity','simulated conductivity','theorical equation')
ylabel('horizontal')
xlabel('m')



%% Error
figure;
for i_sim=1:mm*(kk-1)
    subplot(kk-1,mm,i_sim); pcolor(grid{end}.x,grid{end}.y, Y{end}.m{i_sim}-X_true); 
    title(['Error in Simmulated ',parm.unit, ' Y - X_true']); xlabel('x[m]');ylabel('y[m]');
    shading flat;colorbar;axis tight;set(gca,'Ydir','reverse');
end


%% SuperBlockGrid
i=round(k.sb.nx/3); j=round(k.sb.ny/2);
windows=false(k.sb.ny,k.sb.nx);
for u=1:length(k.el_X_s) %.. look at all point...
    for q=1:4 %... and assign it to the corresponding quadrant
        if i+k.qs(q,1)*k.el_X_s(u)<=k.sb.nx && j+k.qs(q,2)*k.el_Y_s(u)<=k.sb.ny && i+k.qs(q,1)*k.el_X_s(u)>=1 && j+k.qs(q,2)*k.el_Y_s(u)>=1% check to be inside the grid
            windows(i+k.qs(q,1)*k.el_X_s(u), j+k.qs(q,2)*k.el_Y_s(u))=true;
        end
    end
end

figure; hold on
imagesc(k.sb.x,k.sb.y,windows) % windows (yellow)
mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.nx+1,k.sb.ny+1),'EdgeColor','k','facecolor','none') % mesh
plot(X.x, X.y,'d')
plot(k.sb.x(j), k.sb.y(i),'or')
plot(X.x(k.sb.mask(i,j,:)), X.y(k.sb.mask(i,j,:)),'x')
plot([X.x'; k.sb.x(X.sb_x)],[X.y'; k.sb.y(X.sb_y)])
axis tight;
xlabel('x [m]');ylabel('y [m]')
legend('Super Grid','Hard data',['Center of grid ' num2str(i) ';' num2str(j)],['Hard data selected with grid ' num2str(i) ';' num2str(j)])



%% Kernel Estimation
Zdx=nan(size(X.d));
for i=1:X.n
   Zdx(i)=Z.d(X.y(i)==Z.y, X.x(i)==Z.x);
end
figure; hold on;
imagesc(kernel.x, kernel.y, kernel.dens)
plot(Zdx, X.d, 'xk');
xlabel('Primary variable X'); xlabel('Secondary variable Y');
axis tight; colorbar



%% Nscore
figure;
x=195;
h1=subplot(1,2,1);hold on
ecdf(X.d);xlabel('x-original');ylabel('CDF(x)'); xlim([kernel.y(1) kernel.y(end)])
plot(x,Nscore.T_F(x),'.r', 'MarkerSize',40)
line([x x],[0 Nscore.T_F(x)],'Color','k')
line([x kernel.y(end)],[Nscore.T_F(x) Nscore.T_F(x)],'Color','k')
ax = gca; ax.XTick = sort(unique([ax.XTick ,x])); h1.Position(3)=h1.Position(3)+0.03;

h2=subplot(1,2,2);hold on
plot(-5:0.1:5,normcdf(-5:0.1:5))
plot(Nscore.forward(x), Nscore.T_F(x),'.r', 'MarkerSize',40)
line([-5 Nscore.forward(x)],[Nscore.T_F(x) Nscore.T_F(x)],'Color','k')
line([Nscore.forward(x) Nscore.forward(x)],[0 Nscore.T_F(x)],'Color','k')
xlabel('x-Normal Score Transform');ylabel('Standard Normal CDF(x)'); xlim([-5 5])
ax = gca; ax.XTick = sort([ax.XTick ,Nscore.forward(x)]);ax.YAxisLocation='right';
h2.Position(1)=h2.Position(1)-0.03;

% Figure 2
mu= -6;
sigma=0.3;
idx=1:10:kernel.n;

figure;
subplot(2,2,1);hold on
ecdf(X.d);xlabel('x');ylabel('cdf(x)'); xlim([kernel.y(1) kernel.y(end)])
plot(kernel.y(idx),Nscore.T_F(kernel.y(idx)),'.r', 'MarkerSize',20)
for i=idx
    line([kernel.y(i) kernel.y(i)],[0 Nscore.T_F(kernel.y(i))],'Color',[0.4,0.4,0.4])
    line([kernel.y(i) kernel.y(end)],[Nscore.T_F(kernel.y(i)) Nscore.T_F(kernel.y(i))],'Color',[0.4,0.4,0.4])
end

subplot(2,2,2);hold on
plot(-5:0.1:5,normcdf(-5:0.1:5))
xlabel('x');ylabel('Std Normal cdf(x)'); xlim([-5 5])
for i=idx
    line([-5 norminv(Nscore.T_F(kernel.y(i)))],[Nscore.T_F(kernel.y(i)) Nscore.T_F(kernel.y(i))],'Color',[0.4,0.4,0.4])
    line([norminv(Nscore.T_F(kernel.y(i))) norminv(Nscore.T_F(kernel.y(i)))],[0 Nscore.T_F(kernel.y(i))],'Color',[0.4,0.4,0.4])
end
plot(norminv(Nscore.T_F(kernel.y(idx))),Nscore.T_F(kernel.y(idx)),'.r', 'MarkerSize',20)

subplot(2,2,3);hold on
hist(X.d);xlabel('x');ylabel('pdf(x)'); xlim([kernel.y(1) kernel.y(end)])
ax = gca;ax.YAxisLocation='right';set(gca,'Ydir','reverse'); 

subplot(2,2,4); hold on;
plot(-5:0.1:5,normpdf(-5:0.1:5,mu,sigma));
plot(Nscore.forward(kernel.y(idx)),normpdf(Nscore.forward(kernel.y(idx)),mu,sigma),'.r', 'MarkerSize',20)
for i=idx
    line([Nscore.forward(kernel.y(i)) Nscore.forward(kernel.y(i))],[0 normpdf(Nscore.forward(kernel.y(i)),mu,sigma)],'Color','k')
end

for i=1:length(idx)-1
    u=linspace(Nscore.forward(kernel.y(idx(i))),Nscore.forward(kernel.y(idx(i+1))),20);
    area([u  Nscore.forward(kernel.y(idx(i+1))) Nscore.forward(kernel.y(idx(i))) ],...
        [normpdf(u,mu,sigma)  0 0])
end
set(gca,'Ydir','reverse'); xlim([-5 5])















