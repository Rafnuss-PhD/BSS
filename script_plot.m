
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
imagesc(grid_gen.x, grid_gen.y, sigma_true); 
title('True Porosity \rho_{true}'); 
xlabel('x[m]'); ylabel('y [m]'); colorbar;set(gca,'Ydir','reverse');
subplot(4,1,2); hold on; 
imagesc(grid_gen.x, grid_gen.y, log10(K_true)); plot(K.x, K.y, 'or'); 
title('Log True Hydraulic Conudctivity K_{true} and sampled point location K'); 
xlabel('x[m]'); ylabel('y [m]'); colorbar; set(gca,'Ydir','reverse'); axis tight equal
subplot(4,1,3);
imagesc(grid_gen.x, grid_gen.y, sigma_true); hold on; plot(sigma.x, sigma.y, 'or'); 
title('True Electrical Conductivity \sigma_{true} and sampled point location g');
xlabel('x[m]'); ylabel('y [m]');  colorbar; set(gca,'Ydir','reverse'); caxis(caxis_lim);axis equal
subplot(4,1,4); 
imagesc(grid_gen.x, grid_gen.y, Sigma.d); 
title('Inverted Electrical Conductivity \Sigma_{true}'); 
xlabel('x[m]'); ylabel('y [m]');colorbar; set(gca,'Ydir','reverse'); caxis(caxis_lim); axis equal


figure;
subplot(3,1,1); 
imagesc(grid_gen.x,grid_gen.y,sigma_true);hold on; plot(sigma.x, sigma.y, 'or'); 
xlabel('x[m]'); ylabel('y [m]'); colorbar;set(gca,'Ydir','reverse'); caxis(caxis_lim);
title('True Electrical Conductivity [mS/m] and well location');
subplot(3,1,2); 
imagesc(grid_gen.x,grid_gen.y,Sigma.d); 
xlabel('x[m]'); ylabel('y [m]');  colorbar;set(gca,'Ydir','reverse'); caxis(caxis_lim);
title('Electrical Conductivity Tomography [mS/m]');
subplot(3,1,3); 
pcolor(grid_gen.x,grid_gen.y,Sigma.std); shading flat; 
xlabel('x[m]'); ylabel('y [m]'); colorbar; set(gca,'Ydir','reverse');
title('Electrical Conductivity Tomography standard deviation [mS/m]'); 
% subplot(4,1,4); 
% imagesc(grid_gen.x,grid_gen.y,sigma_true-Sigma.d); 
% xlabel('x[m]'); ylabel('y [m]'); title('rho_{true}-G');colorbar;set(gca,'Ydir','reverse');


figure;
subplot(2,1,1); 
imagesc(grid_gen.x,grid_gen.y,log(K_true));hold on; plot(sigma.x, sigma.y, 'or'); 
xlabel('x[m]'); ylabel('y [m]'); colorbar;set(gca,'Ydir','reverse'); 
title('True hydraulic conductivity [Log m/s]'); axis equal tight
caxis([kern.axis_prim(1) kern.axis_prim(end)])
subplot(2,1,2); 
imagesc(grid_gen.x,grid_gen.y,Sigma.d); 
xlabel('x[m]'); ylabel('y [m]');  colorbar;set(gca,'Ydir','reverse');
title('Electrical conductivity from ERT [mS/m]'); axis equal tight
caxis([kern.axis_sec(1) kern.axis_sec(end)])

% Histogram
figure;
nbins_all = grid_gen.nxy/1000;
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




%% BSGS
% Multi-grid
addpath(genpath('./.'))
sub_n = min([parm.n_realisation,5])+1;
if sub_n<4, kk=sub_n;mm=1;
else kk=ceil(sub_n/2); mm=2; 
end

grid=grid_gen;
X_true = log(K_true);
Y = Res;
X = Klog;
Z = Sigma;
kernel = kern;
kernel.y = kern.axis_prim;
kernel.x = kern.axis_sec;


figure(2);clf;
caxis_limm = [min([X_true(:);Y.m{end}(:)]) max([X_true(:);Y.m{end}(:)])];

subplot(kk,mm,1); hold on; 
pcolor(grid_gen.x,grid_gen.y,X_true);
plot(X.x,X.y,'or')
shading flat;colorbar; caxis(caxis_limm); axis tight;set(gca,'Ydir','reverse'); xlabel('x[m]'); ylabel('y[m]');
title([ 'True ',parm.unit ' X_{true}: \mu=' num2str(mean2(X_true)) ' | \sigma='   num2str(std2(X_true)) ' and X_{sampled}: \mu=' num2str(mean(X.d)) ' | \sigma='   num2str(std(X.d))]);
for i_sim=2:mm*kk
    subplot(kk,mm,i_sim); 
    pcolor(Y.x,Y.y,Y.m{i_sim-1+9});
    shading flat;colorbar;caxis(caxis_limm);axis tight;set(gca,'Ydir','reverse'); xlabel('x[m]');ylabel('y[m]');
    title(['Simmulated ', parm.unit, ': \mu=' num2str(mean2(Y.m{i_sim-1})) ' | \sigma=' num2str(std2(Y.m{i_sim-1})) ]);
end
if isfield(parm, 'savefig') && parm.savefig
    filename=['result/', parm.familyname, 'simulations_', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.fig'];
    savefig(filename)
end

% Histogramm
bins=linspace(min(min([Y.m{:}])), max(max([Y.m{:}])), 25);
ff=zeros(numel(bins),parm.n_realisation);
for i_sim=1:parm.n_realisation
    f=hist(Y.m{i_sim}(:),bins); 
    ff(:,i_sim) = f/trapz(bins,f);
end
figure(3);clf;  hold on; 
plot(bins,mean(ff,2),'--k')
plot(bins,max(ff,[],2),'k')
plot(bins,min(ff,[],2),'k')
h1=fill( [bins fliplr(bins)],  [max(ff,[],2)' fliplr(min(ff,[],2)')], [.4 .4 .4],'DisplayName','min and max range '); alpha(.25);
% stairs(x,f/trapz(x,f),'DisplayName',sprintf('Simulation %d | \\mu=%.2f, \\sigma=%.2f',i_sim,mean(Y.m{i_sim}(:)), std(Y.m{i_sim}(:))));
f=hist(X_true(:),bins); h2=plot(bins,f/trapz(bins,f),'linewidth',2,'DisplayName',sprintf('True : \\mu=%.1f | \\sigma=%.1f',mean(X_true(:)), std(X_true(:))));
[f,x]=hist(X.d(:),bins); h3=plot(x,f/trapz(x,f),'linewidth',2,'DisplayName',sprintf('Sampled| \\mu=%.1f | \\sigma=%.1f',mean(X.d(:)), std(X.d(:))));

legend([h2, h3 h1]); axis tight
ylabel('Histogram');xlabel('Electrical conductivity [mS/m]');


if isfield(parm, 'savefig') && parm.savefig
    filename=['result/', parm.familyname, 'histo_', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.fig'];
    savefig(filename)
end


%% Variogram
[gamma_x_s, gamma_y_s] = variogram_gridded_perso(X_true);
for i_realisation=1:parm.n_realisation
    [gamma_x_y{i_realisation}, gamma_y_y{i_realisation}] = variogram_gridded_perso(Y.m_ns{i_realisation});
end


figure(4);clf; 
subplot(1,2,1);
hold on
id = grid_gen.x<parm.k.covar(1).range(1)*parm.k.wradius;
for i_realisation=1:parm.n_realisation
    h1=plot(grid_gen.x(id),gamma_x_y{i_realisation}(id),'Color', [.4 .4 .4]);
end
h3=plot(grid_gen.x(id),gamma_x_s(id)./std(X_true(:)),'k','linewidth',2);
h2=plot(grid_gen.x(id),(1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1))),'--k','linewidth',3);
legend([h2, h3, h1],'Theorical model C(h)','True','simulation(s)');
ylabel('Horizontal'); xlabel('m');axis tight

% filename=['result/1-Generale/vario_h_', datestr(now,'yyyy-mm-dd_HH-MM-SS')];
% savefig(filename);set(gcf, 'Color', 'w');
% export_fig([filename '.eps'])

% 
subplot(1,2,2);
hold on
id= grid_gen.y<parm.k.covar(1).range(2)*1.2;

for i_realisation=1:parm.n_realisation
    h1=plot(grid_gen.y(id),gamma_y_y{i_realisation}(id),'Color', [.4 .4 .4]);
end
h3=plot(grid_gen.y(id),gamma_y_s(id),'k','linewidth',3);
h2=plot(grid_gen.y(id),(1-parm.k.covar(1).g(grid_gen.y(id)/parm.k.covar(1).range(2))),'--k','linewidth',3);
legend([h2, h3, h1],'Theorical model C(h)','True','simulation(s)','Location','northwest')
ylabel('Vertical'); xlabel('m'); axis tight



%% Several field statistic
YY=nan(parm.n_realisation,numel(Y.m{1}));
for i_realisation=1:parm.n_realisation
    YY(i_realisation,:) = Y.m{i_realisation}(:);
end

YY_mean=reshape(mean(YY),size(Y.m{1},1),size(Y.m{1},2));
YY_std=reshape(std(YY),size(Y.m{1},1),size(Y.m{1},2));

figure(5);clf; hold on;
subplot(3,1,1); imagesc(grid.x,grid.y,X_true); colorbar; title('true primary'); c_axis=caxis;
subplot(3,1,2); imagesc(grid.x,grid.y,YY_mean); colorbar;title('Mean of realisations'); caxis(c_axis);
subplot(3,1,3); imagesc(grid.x,grid.y,YY_std); colorbar; title('Std of realisations')



%% Error
figure(6);clf
for i_sim=1:mm*(kk-1)
    subplot(kk-1,mm,i_sim); pcolor(grid.x,grid.y, Y.m{i_sim}-X_true); 
    title(['Error in Simmulated ',parm.unit, ' Y - X_{true}']); xlabel('x[m]');ylabel('y[m]');
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
mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.ny+1, k.sb.nx+1),'EdgeColor','k','facecolor','none') % mesh
plot(X.x, X.y,'d')
plot(k.sb.x(j), k.sb.y(i),'or')
plot(X.x(k.sb.mask{i,j}), X.y(k.sb.mask{i,j}),'x')
% plot([X.x'; k.sb.x(X.sb_x)],[X.y'; k.sb.y(X.sb_y)])
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
ylabel('Primary variable X'); xlabel('Secondary variable Y');
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




%% 










