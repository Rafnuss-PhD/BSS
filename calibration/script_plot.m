
%% 
figure('Position', [100, 100, 700, 500]);
K_true_log=log(K_true);

subplot(2,3,1); hold on; 
[X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
jpdf = ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]);
Resjpdf=reshape(mean(E2,2)+jpdf,size(X,1),size(X,2));
imagesc(kern.axis_sec, kern.axis_prim,Resjpdf);
Zdx=nan(size(Klog.d));
for i=1:Klog.n
   Zdx(i)=Sigma.d(Klog.y(i)==Sigma.y, Klog.x(i)==Sigma.x);
end
plot(Zdx, Klog.d, 'xk');
axis tight; %colorbar; 
title('(a) Joint PDF');
xlabel('\sigma')
ylabel('K')


subplot(2,3,2);  hold on; 
id = grid_gen.x<parm.k.covar(1).range(1)*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)))';
a = reshape(Nscore.forward(K_true_log(:)),grid_gen.ny,grid_gen.nx);
gamma_x_s = variogram_gridded_perso( a./std(a(:)) );
h1=plot(grid_gen.x(id),E1+repmat(Gamma_t,1,parm.n_realisation),'Color', [.4 .4 .4]);
h3=plot(grid_gen.x(id),gamma_x_s(id),'k','linewidth',2);
h2=plot(grid_gen.x(id),Gamma_t,'--k','linewidth',3);
axis tight; caxis([0 .3]);
xlabel('h'); ylabel('C(h)')
title('Covariance function')

subplot(2,3,3)
YY=nan(parm.n_realisation,numel(Res.m{1}));
for i_realisation=1:parm.n_realisation
    YY(i_realisation,:) = Res.m{i_realisation}(:)-K_true_log(:);
end
% YY_mean=sqrt(reshape(mean(YY),Res.ny,Res.nx).^2);
% imagesc(grid_gen.x,grid_gen.y,YY_mean);
% xlabel('x'); ylabel('y'); axis equal tight;
histogram(YY(YY(:)~= 0)); axis tight
title('Histogram of error')

subplot(2,3,[4 6])
imagesc(grid_gen.x,grid_gen.y,Res.m{1});
caxis([kern.axis_prim(1) kern.axis_prim(end)])
xlabel('x'); ylabel('y'); axis equal tight; %  colorbar;

set(gca, 'Color', 'none');
% saveas(gcf, 'result-BSS/figure/Res11.png');
% saveas(gcf, 'result-BSS/figure/Res11.fig');
% saveas(gcf, 'result-BSS/figure/Res11.eps');

% export_fig result-BSS/figure/colorbar  -png -m2 -transparent

%%

files={'Res0','Res1','Res11','Res5'};
names={'w_X=1 and w_Z=0\nkriging of K', 'w_X=0 and w_Z=1\nsequential sim.','w_X=w_Z=1$\ntraditional BSS','w_X=w_Z=.5\nnew BSS'};
n=numel(files);

Zdx=nan(size(Klog.d));
for i=1:Klog.n
    Zdx(i)=Sigma.d(Klog.y(i)==Sigma.y, Klog.x(i)==Sigma.x);
end
[X,Y] = meshgrid(kern.axis_sec, kern.axis_prim);
jpdf = ksdensity([sigma.d(:) Klog.d],[X(:),Y(:)]);
id = grid_gen.x<parm.k.covar(1).range(1)*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x(id)/parm.k.covar(1).range(1)))';
a = reshape(Nscore.forward(K_true_log(:)),grid_gen.ny,grid_gen.nx);
gamma_x_s = variogram_gridded_perso( a./std(a(:)) );

figure;% ('Position', [100, 100, 700, 500]);  
for u=1:n
    load(['Y:\BSGS\result-BSS\',files{u},'.mat'], 'E1', 'E2', 'parm')
    subplot(2,n,u); hold on;
    imagesc(kern.axis_sec, kern.axis_prim,reshape(mean(E2,2)+jpdf,size(X,1),size(X,2)));
    scatter(Zdx, Klog.d, '.w');
    axis tight; %colorbar;
    title(names{u});
    xlabel('\sigma [mS/m]')
    ylabel('K [log m/s]')
    
    subplot(2,n,n+u);  hold on;
    h1=plot(grid_gen.x(id),E1+repmat(Gamma_t,1,parm.n_realisation),'Color', [.4 .4 .4]);
    h3=plot(grid_gen.x(id),gamma_x_s(id),'k','linewidth',2);
    h2=plot(grid_gen.x(id),Gamma_t,'--k','linewidth',3);
    axis tight; ylim([0 1.4]);
    ylabel('C(h)'); xlabel('h [m]')
end

figure; hold on;
for u=1:n
    load(['Y:\BSGS\result-BSS\',files{u},'.mat'], 'Res')
    YY=nan(parm.n_realisation,numel(Res.m{1}));
    for i_realisation=1:parm.n_realisation
        YY(i_realisation,:) = Res.m{i_realisation}(:)-K_true_log(:);
    end
    % YY_mean=sqrt(reshape(mean(YY),Res.ny,Res.nx).^2);
    % imagesc(grid_gen.x,grid_gen.y,YY_mean);
    % xlabel('x'); ylabel('y'); axis equal tight;
    ksdensity(YY(:)); axis tight
    title('Histogram of error')
end
legend(names)

figure; hold on;
for u=1:n
    load(['Y:\BSGS\result-BSS\',files{u},'.mat'], 'Res')
    subplot(n+1,1,u)
    imagesc(grid_gen.x,grid_gen.y,Res.m{1});
    caxis([kern.axis_prim(1) kern.axis_prim(end)])
    xlabel('x'); ylabel('y'); axis equal tight;
    title(names{u})
end
subplot(n+1,1,n+1)
imagesc(grid_gen.x,grid_gen.y,K_true_log)
caxis([kern.axis_prim(1) kern.axis_prim(end)])
xlabel('x'); ylabel('y'); axis equal tight;
colorbar



%% 
figure; 
subplot(4,1,[1 3]); hold on;
files={'ResAB','ResAB3','ResAB2','ResAB4','Res11'};
dx=.006; dy=.0001;
for i=1:numel(files)
    load(['Y:\BSGS\result-BSS\',files{i},'.mat'], 'OF1', 'OF2', 'parm')
    scatter(OF1,OF2)
    xlabel('OF1 (Variogram)')
    ylabel('OF2 (Joint PDF)')
    text(OF1+dx,OF2+dy,strread(num2str(parm.aggr.a),'%s'),'HorizontalAlignment','left');
end
legend({'b=Inf','b=100','b=10','b=5','Initial BSS'})

subplot(4,1,4); hold on;
i_pt_temp=1:grid_gen.nxy;
plot(i_pt_temp,parm.aggr.fx(.5,realmax,i_pt_temp./grid_gen.nxy),'linewidth',2);
plot(i_pt_temp,parm.aggr.fx(.1,100,i_pt_temp./grid_gen.nxy),'linewidth',2);
plot(i_pt_temp,parm.aggr.fx(.15,10,i_pt_temp./grid_gen.nxy),'linewidth',2);
plot(i_pt_temp,parm.aggr.fx(0,5,i_pt_temp./grid_gen.nxy),'linewidth',2);
plot([0 i_pt_temp(end)],[1 1],'linewidth',2);
legend({'a=.5 ; b=Inf', 'a=.1 ; b=100', 'a=.15 ; b=10', 'a=0 ; b=5', 'Initial BSS'})
