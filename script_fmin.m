%% Preparation
load('result-BSS/GEN-Run_1_2017-05-07_14-37');
addpath('/home/rnussbaumer/MATLAB/')
addpath(genpath('./.'));
Nscore = nscore(kern, struct('nscore', 1), 0); %Prim.d, kern.axis_prim, 'pchip', 'pchip', parm.plot.ns
% [~,s_id]=min(bsxfun(@minus,kern.axis_sec,Sigma.d(:)).^2,[],2);
sec_pdf = kern.dens(:,s_id);
sec.pdf = bsxfun(@times, sec_pdf, 1./sum(sec_pdf));
sec.axis = Nscore.forward(kern.axis_prim);

parm.k.covar = gen.covar;
parm.k.covar.range0 = fliplr(gen.covar.range0) ./ [grid_gen.dy grid_gen.dx];
parm.k.covar = kriginginitiaite(parm.k.covar);

parm.k.wradius = 2;
parm.k.lookup = false;
parm.k.nb = 100;
parm.path_random     = 1;
parm.path            = 'linear';
parm.mg=0;

parm.seed_path      = 'shuffle';
parm.seed_search    = 'shuffle';

parm.aggr.method='linear';
parm.n_real=1;

% use the log of hyd. cond.
hd = sampling_pt(struct('x',1:grid_gen.nx,'y',1:grid_gen.ny),log(K_true),1,0);
hd.d = Nscore.forward(hd.d);
hd.scale=nan(hd.n,1);
hd.x=hd.x(:); hd.y=hd.y(:); hd.d=hd.d(:); 

f0=kern.prior ./ sum(kern.prior);
nx = grid_gen.nx;
ny = grid_gen.ny;

k = parm.k;

% BSS function

% 1. Creation of the grid an path
[Y, X] = ndgrid(1:ny,1:nx);


% 2. Define Path
Path = nan(ny,nx);
Path(hd.id) = 0;

rng(parm.seed_path);
if parm.mg
   sx = 1:ceil(log(nx+1)/log(2));
   sy = 1:ceil(log(ny+1)/log(2));
   sn = max([numel(sy), numel(sx)]);
   nb = nan(sn,1);
   start = zeros(sn+1,1);
   dx = nan(sn,1); dy = nan(sn,1);
   path = nan(sum(isnan(Path(:))),1);
   for i_scale = 1:sn
       dx(i_scale) = 2^(numel(sx)-sx(min(i_scale,end)));
       dy(i_scale) = 2^(numel(sy)-sy(min(i_scale,end)));
       [Y_s,X_s] = ndgrid(1:dy(i_scale):ny,1:dx(i_scale):nx); % matrix coordinate
       id = find(isnan(Path(:)) & ismember([Y(:) X(:)], [Y_s(:) X_s(:)], 'rows'));
       nb(i_scale) = numel(id);
       start(i_scale+1) = start(i_scale)+nb(i_scale);
       path( start(i_scale)+(1:nb(i_scale)) ) = id(randperm(nb(i_scale)));
       Path(path( start(i_scale)+(1:nb(i_scale)) )) = start(i_scale)+(1:nb(i_scale));
       % Find the scaloe of hard data.
       hd.scale( ismember([hd.y hd.x], [Y_s(:) X_s(:)],'rows') & isnan(hd.scale)) =i_scale;
   end
else
   id=find(isnan(Path));
   path = id(randperm(numel(id)));
   Path(path) = 1:numel(id);
   dx=1; dy=1; nb = numel(id); start=[0 nb]; sn=1;
end


% 3. Initialization Spiral Search
% Initialize spiral search stuff which don't change
x = ceil( min(k.covar(1).range0(2)*k.wradius, nx));
y = ceil( min(k.covar(1).range0(1)*k.wradius, ny));
[ss_Y, ss_X] = ndgrid(-y:y, -x:x);% grid{i_scale} of searching windows
ss_dist = sqrt( (ss_X/k.covar(1).range(2)).^2 + (ss_Y/k.covar(1).range(1)).^2); % find distence
ss_id_1 = find(ss_dist <= k.wradius); % filter node behind radius.
rng(parm.seed_search);
ss_id_1 = ss_id_1(randperm(numel(ss_id_1)));
[ss_dist_s, ss_id_2] = sort(ss_dist(ss_id_1)); % sort according distence.
ss_X_s=ss_X(ss_id_1(ss_id_2)); % sort the axis
ss_Y_s=ss_Y(ss_id_1(ss_id_2));
ss_n=numel(ss_X_s); %number of possible neigh

if parm.mg
    ss_scale=sn*ones(size(ss_X));
    for i_scale = sn-1:-1:1
        x_s = [-fliplr(dx(i_scale):dx(i_scale):x(end)) 0 dx(i_scale):dx(i_scale):x(end)]+(x+1);
        y_s = [-fliplr(dy(i_scale):dy(i_scale):y(end)) 0 dy(i_scale):dy(i_scale):y(end)]+(y+1);
        ss_scale(y_s,x_s)=i_scale;
    end
    ss_scale_s = ss_scale(ss_id_1(ss_id_2));
else
    ss_scale_s = sn*ones(size(ss_id_2));
end


% 4. Simulation
tik.weight = tic;
NEIGH = nan(nx*ny,k.nb);
% NEIGH_1 = nan(nx*ny,k.nb);
% NEIGH_2 = nan(nx*ny,k.nb);
LAMBDA = nan(nx*ny,k.nb);
S = nan(nx*ny,1);

k_nb = k.nb;
k_covar_c0 = sum([k.covar.c0]);
XY_i=[Y(path) X(path)];

for i_scale = 1:sn
    ss_id = find(ss_scale_s<=i_scale);
    ss_XY_s_s = [ss_Y_s(ss_id) ss_X_s(ss_id)];
    ss_dist_s_s = ss_dist_s(ss_id);

    
    % Remove hard data which are on the current scale
    hd_XY_s = [hd.y(hd.scale>i_scale) hd.x(hd.scale>i_scale)];
    
    for i_pt = start(i_scale)+(1:nb(i_scale))
        
        % Compute distance to the hard data
        hd_XY_d = bsxfun(@minus,hd_XY_s,XY_i(i_pt,:));
        hd_XY_d = hd_XY_d(hd_XY_d(:,1)<k.covar(1).range(1)*k.wradius &  hd_XY_d(:,2)<k.covar(1).range(2)*k.wradius,:);
        hd_dist=zeros(size(hd_XY_d,1),1);
        for i=1:numel(k.covar)
            hd_dist=hd_dist+sqrt(sum((hd_XY_d*k.covar(i).cx).^2,2));
        end
        
        [~, ss_hd_id] = sort( [ hd_dist; ss_dist_s_s]);
        tmp = [hd_XY_d; ss_XY_s_s];
        ss_hd_XY_s_s = tmp(ss_hd_id,:);
        
        % Neighborhood search
        n=0;
        neigh=nan(k_nb,1);
        NEIGH_1 = nan(k.nb,1);
        NEIGH_2 = nan(k.nb,1);
        for nn = 2:size(ss_hd_XY_s_s,1) % 1 is the point itself... therefore unknown
            ijt = XY_i(i_pt,:) + ss_hd_XY_s_s(nn,:);
            if ijt(1)>0 && ijt(2)>0 && ijt(1)<=ny && ijt(2)<=nx
                if Path(ijt(1),ijt(2)) < i_pt % check if it,jt exist
                    n=n+1;
                    neigh(n) = nn;
                    NEIGH_1(n) = ijt(1);
                    NEIGH_2(n) = ijt(2);
                    if n >= k_nb
                        break;
                    end
                end
            end
        end
        
        
        % Covariance computation
        if n==0
            S(i_pt) = k_covar_c0;
        else
            NEIGH(i_pt,:) = NEIGH_1 + (NEIGH_2-1)* ny;
            D = pdist([0 0; ss_hd_XY_s_s(neigh(1:n),:)]*k.covar.cx);
            C = k.covar.g(D);
            
            if n==1
                a0_C = C;
                ab_C = 1;
            else
                a0_C = C(1:n)';
                % Equivalent to : squareform(C(n+1:end));
                ab_C = diag(ones(n,1))*0.5;
                ab_C(tril(true(n),-1)) = C(n+1:end);
                ab_C = ab_C + ab_C';
            end

            % Weights and variance error
            l = ab_C \ a0_C;
            LAMBDA(i_pt,:) = [l; nan(k.nb-n,1)];
            S(i_pt) = k_covar_c0 - l'*a0_C;
        end
    end
    % disp(['scale: ' num2str(i_scale) '/' num2str(sn)])
end
t.weight = toc(tik.weight);
disp(['Weights Computed in ' num2str(t.weight*60)] )


% Prepare OF
clear id;
id.x = (1:ceil(k.covar(1).range0(2)*k.wradius))';
id.y = (1:ceil(k.covar(1).range0(1)*k.wradius))';
Gamma_id.x = (1-k.covar.g([0 ; id.x(1:end-1)]*k.covar.cx(4)))';
Gamma_id.y = (1-k.covar.g([0 ; id.y(1:end-1)]*k.covar.cx(1)))';

XY = kern.XY;
XY(:,2)= Nscore.forward(XY(:,2));
Sigma_d = Sigma.d(:);
dens = kern.dens(:)./sum(kern.dens(:));

inve = @(x)  reshape(Nscore.inverse(reshape(x,numel(x),1)),size(x,1),size(x,2));

% Evaluation of the true
KlogNscore = reshape(Nscore.forward(log(K_true(:))), grid_gen.ny,  grid_gen.nx);
r = KlogNscore;
% r = (r-mean(r(:)))/std(r(:));
[gamma_x_s, gamma_y_s] = variogram_gridded_perso(r);
dens_r =ksdensity([Sigma_d r(:)],XY);
E2t=dens_r./sum(dens_r) - dens;
OF1t = sqrt(sum((((1-Gamma_id.x(:)).*(gamma_x_s(id.x)-Gamma_id.x(:))).^2))) + sqrt(sum((((1-Gamma_id.y(:)).*(gamma_y_s(id.y)-Gamma_id.y(:))).^2))) ;
OF2t = sqrt(mean(E2t.^2));


%% Generated Data
figure(1);
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






%% Cst Weight

T0 = [1 0;0 1;.5 .5; 1 1; .9 .6; 1.1 .3];
parm.n_real=48*size(T0,1)*2;
parm.aggr.method='cst';
[Res, out, OF1, OF2, E1_x, E1_y, E2] = fmin(T0,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
% Changhe the sum or the prior
T02=[1 1];
parm2= parm; parm2.n_real=48*size(T02,1)*2;
f0t = f0*0+1;
[Res2, out2, OF12, OF22, E1_x2, E1_y2, E22] = fmin(T02,parm2,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0t,id,kern,Gamma_id,Sigma_d,XY,dens,hd);



figure(2); % Exemple of each
subplot(size(T0,1)+2,1,1);
imagesc(grid_gen.x,grid_gen.y,log(K_true));
caxis([kern.axis_prim(1) kern.axis_prim(end)])
subplot(size(T0,1)+2,1,2);
imagesc(grid_gen.x,grid_gen.y,inve(Res2(:,:,1)));
caxis([kern.axis_prim(1) kern.axis_prim(end)])
for i_res = 1:size(T0,1)
    subplot(size(T0,1)+2,1,i_res+2);
    imagesc(grid_gen.x,grid_gen.y, inve(Res(:,:,i_res)))
    caxis([kern.axis_prim(1) kern.axis_prim(end)])
end
colorbar;

figure(3); % Joint PDF
subplot(2,3,1); hold on;
%imagesc(kern.axis_sec, kern.axis_prim, kern.dens); 
imagesc(kern.axis_sec, kern.axis_prim, reshape(E2t+dens,numel(kern.axis_prim),numel(kern.axis_sec))); 
xlim([2 12]); ylim([-10 -4]);axis square;caxis([0 15]*10^-4)
subplot(2,3,2); hold on;
imagesc(kern.axis_sec, kern.axis_prim, reshape(mean(E22,2)+dens,numel(kern.axis_prim),numel(kern.axis_sec)));
xlim([2 12]); ylim([-10 -4]);axis square;caxis([0 15]*10^-4)
for i_res = 1:size(T0,1)
    subplot(2,3,i_res+2); hold on;
    ii_t=i_res:size(T0,1):parm.n_real;
    imagesc(kern.axis_sec, kern.axis_prim, reshape(mean(E2(:,ii_t),2)+dens,numel(kern.axis_prim),numel(kern.axis_sec)));
    xlim([2 12]); ylim([-10 -4]); axis square;  caxis([0 15]*10^-4)
end


figure(4); % Variogram
col=get(groot,'DefaultAxesColorOrder');
subplot(2,1,1); hold on;
h(i_res+1)=shadedErrorBar(grid_gen.x(id.x),mean(E1_x2,2)+Gamma_id.x',[max(E1_x2,[],2)-mean(E1_x2,2)  abs(min(E1_x2,[],2)-mean(E1_x2,2))],{'color',col(1,:)},1);
for i_res = 1:size(T0,1)
    ii_t=i_res:size(T0,1):parm.n_real;
%     for i_realisation=1:parm.n_realisation
%         % gamma_x = variogram_gridded_perso(Res.m_ns{i_realisation});
%         plot(grid_gen.x(id),E1(id, i_realisation)+Gamma_t_id,'Color', [.4 .4 .4]);
%     end
    h(i_res) = shadedErrorBar(grid_gen.x(id.x),mean(E1_x(:,ii_t),2)+Gamma_id.x',[max(E1_x(:,ii_t),[],2)-mean(E1_x(:,ii_t),2)  abs(min(E1_x(:,ii_t),[],2)-mean(E1_x(:,ii_t),2))],{'color',col(i_res+1,:)},1);
end
h2(1)=plot(grid_gen.x(id.x),gamma_x_s(id.x),'k','linewidth',2);
h2(2)=plot(grid_gen.x(id.x),Gamma_id.x,'--k','linewidth',3);
ylabel('C(h)'); xlabel('h_x [m]'); axis tight
subplot(2,1,2); hold on;
h(i_res+1)=shadedErrorBar(grid_gen.y(id.y),mean(E1_y2,2)+Gamma_id.y',[max(E1_y2,[],2)-mean(E1_y2,2)  abs(min(E1_y2,[],2)-mean(E1_y2,2))],{'color',col(1,:)},1);
for i_res = 1:size(T0,1)
    ii_t=i_res:size(T0,1):parm.n_real;
%     for i_realisation=1:parm.n_realisation
%         % gamma_y = variogram_gridded_perso(Res.m_ns{i_realisation});
%         plot(grid_gen.y(id),E1(id, i_realisation)+Gamma_t_id,'Color', [.4 .4 .4]);
%     end
    h(i_res) = shadedErrorBar(grid_gen.y(id.y),mean(E1_y(:,ii_t),2)+Gamma_id.y',[max(E1_y(:,ii_t),[],2)-mean(E1_y(:,ii_t),2)  abs(min(E1_y(:,ii_t),[],2)-mean(E1_y(:,ii_t),2))],{'color',col(i_res+1,:)},1);
end
h2(1)=plot(grid_gen.y(id.y),gamma_y_s(id.y),'k','linewidth',2);
h2(2)=plot(grid_gen.y(id.y),Gamma_id.y,'--k','linewidth',3);
ylabel('C(h)'); ylabel('h_y [m]'); axis tight


figure(5); % OF plot
hold on;
scatter(OF12,OF22,'filled')
scatter(OF1,OF2,'filled')
scatter(OF1t,OF2t,'filled')


%% Step vs Variable cst

parm.aggr.method='cst';

[d.T1,d.T2] = meshgrid(0:.05:.5,1);
parm.n_real=numel(d.T1)*48*2;
[~,d.out,d.OF1, d.OF2] = fmin([0.9 0.6],parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);


% OFx,OFy, OF2
figure(11); clf;
subplot(1,3,1);hold on;
imagesc(d.T1(1,:),d.T2(:,1),reshape(d.OF1,numel(d.T1(1,:)),numel(d.T2(:,1))));
axis square tight; caxis(OF1_range)
subplot(1,3,2); hold on;
imagesc(d.T1(1,:),d.T2(:,1),reshape(d.OF2,numel(d.T1(1,:)),numel(d.T2(:,1))));
axis square tight;caxis(OF2_range)
subplot(1,3,3); hold on;
imagesc(d.T1(1,:),d.T2(:,1),reshape(d.out,numel(d.T1(1,:)),numel(d.T2(:,1))));
% [~,min_out] = min(d.out);;plot([d.T1(min_out) d.T1(min_out) 0],[0 d.T2(min_out) d.T2(min_out)],'r','linewidth',2);



OF1_range= [0.0141 1.6420];
OF2_range=1.0e-04 *[0.3398 0.8857];
d.out1 = (d.OF1-OF1_range(1))./(range(OF1_range));
d.out2 = (d.OF2-OF2_range(1))./(range(OF2_range));
w=0:.001:1; clear min_out;
for i_w=1:numel(w)
    out = w(i_w)*d.out1 + (1-w(i_w))*d.out2;
    [~,min_out(i_w)] = min(out);
end
figure(16)
subplot(1,2,1);hold on;
imagesc(d.T1(1,:),d.T2(:,1),reshape(d.OF1,numel(d.T1(1,:)),numel(d.T2(:,1))));
plot(d.T1(min_out), d.T2(min_out),'.-r')
axis square tight; caxis(OF1_range)
subplot(1,2,2); hold on;
imagesc(d.T1(1,:),d.T2(:,1),reshape(d.OF2,numel(d.T1(1,:)),numel(d.T2(:,1))));
plot(d.T1(min_out), d.T2(min_out),'.-r')
axis square tight;caxis(OF2_range)




db.T=[d.T1(unique(min_out)') d.T2(unique(min_out)')];
parm.n_real=size(db.T,1)*48;
[db.Res, db.out, db.OF1, db.OF2, db.E1_x, db.E1_y, db.E2] = fmin(db.T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);



figure(17); % Exemple of each
for i_t=2:size(db.T,1)
    subplot(size(db.T,1)-1,2,(i_t-2)*2+1); 
    imagesc(grid_gen.x,grid_gen.y,inve(db.Res(:,:,i_t)));
    caxis([kern.axis_prim(1) kern.axis_prim(end)])
    subplot(size(db.T,1)-1,2,(i_t-2)*2+2); 
    imagesc(grid_gen.x,grid_gen.y,inve(db.Res(:,:,i_t+size(db.T,1))));
    caxis([kern.axis_prim(1) kern.axis_prim(end)])
end


%% Step
parm.aggr.method='step';

T0=(0.05:.01:.2)'; T0=[T0 ones(numel(T0),1)];
parm.n_real=size(T0,1)*48;
[cst.Res,cst.out,cst.OF1, cst.OF2] = fmin(T0,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);


figure(5); % OF plot
hold on;
scatter(OF12,OF22,'filled')
scatter(OF1,OF2,'filled')
scatter(OF1t,OF2t,'filled')
scatter(db.OF1, db.OF2,'filled')
scatter(cst.OF1, cst.OF2,'filled')


T1=[0:.01:.1 .15:.05:.5 .6:.1:1]'; 
T2=.5:.1:2;
[T1,T2]=meshgrid(T1,T2);
T0=[T1(:),T2(:)];
parm.n_real=size(T0,1)*48;
[cststep.Res,cststep.out,cststep.OF1, cststep.OF2] = fmin(T0,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);

figure(5); % OF plot
hold on;
scatter(OF12,OF22,'filled')
scatter(OF1,OF2,'filled')
scatter(OF1t,OF2t,'filled')
scatter(db.OF1, db.OF2,'filled')
scatter(cst.OF1, cst.OF2,'filled')
scatter(cststep.OF1, cststep.OF2,[],T2(:),'filled')


%% MATLAB CALIBRATION

Run Fminc
parm.n_real=48;

parm.aggr.method='linear';

OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);

OF([2]');


T0 = [.5 .5]; % T_1, T_2 w_1 w_2,
T0 = [0.9    1   0.9 0.1 ]; % T_1, T_2 w_1 w_2,
T0 = [0.3    0.6    .3    .6]; % T_1, T_2 w_min w_max,
out = OF(T0);

out = OF(T0);
out = OF([3 NaN 0 NaN]);




A=[1 -1 0 0 ; 0 0 -1 1];
b=[0;0];
lb=[0 0 0 0];
ub=[1 1 1 1];

options = optimoptions('fmincon','Display','iter-detailed','Diagnostics','on','MaxFunctionEvaluations',100,...
'PlotFcn',{@optimplotx, @optimplotfval , @optimplotstepsize  },'UseParallel',false);
x = fmincon(OF,T0,[],[],[],[],lb,ub,[],options)


options = optimoptions(@particleswarm,'PlotFcn',@pswplotranges,'Display','iter','UseVectorized',true);
x = particleswarm(OF,4,lb,ub);


options = optimoptions(@ga,'PlotFcn',{@gaplotbestf , @gaplotbestindiv, @gaplotexpectation, @gaplotrange},'Display','iter','UseVectorized',true);
x = ga(OF,4,A,b,[],[],lb,ub,[],options);
x = [0.4782    0.8233    0.6332    0.3947];

options = optimoptions(@simulannealbnd,'PlotFcn',{@saplotbestf, @saplotbestx, @saplotf, @saplotx, @saplotstopping , @saplottemperature},'Display','iter','PlotInterval',10);
x = simulannealbnd(OF,x,lb ,ub,options);
x = [0.5556    0.3521    0.5614    0.5666];

options = optimoptions(@patternsearch,'PlotFcn',{@psplotbestf, @psplotmeshsize, @psplotfuncount, @psplotbestx},'Display','iter','PlotInterval',1);
x = patternsearch(OF,T0,A,b,[],[],lb,ub,[],options);


%% eXHAUSITVE SEARCH
n_pool=48;
parpool(n_pool);
parm.n_real=numel(d01.T1)*48;

parm.aggr.method='cst';
parm.n_real=4*48*6;
OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
[d.OF1, d.OF2] = OF([0 1;1 0;.5 .5; 1 1]);
OF1_range= [0.0516 0.2652]; %[d.OF1(2) d.OF1(1)]; %
OF2_range=1.0e-04 *[0.3251    0.7604];  %[d.OF2(1) d.OF2(2)]; % 


parm.aggr.method='cst';
OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
[d.T1,d.T2] = meshgrid(0:.5:2,0:.5:2);
[d.OF1, d.OF2] = OF([d.T1(:) d.T2(:)]);
[d01.T1,d01.T2] = meshgrid(0:.2:1,0:.2:1);
[d01.OF1, d01.OF2] = OF([d01.T1(:) d01.T2(:)]);
[d68.T1,d68.T2] = meshgrid(0.6:.05:1,0.2:.05:.6);
[d68.OF1, d68.OF2] = OF([d68.T1(:) d68.T2(:)]);

save('cst2WeightsVariable','parm','d','d01','d02','d68');

parm.aggr.method='step';
OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
[d2.OF1, d2.OF2, d2.OF1std, d2.OF2std ] = OF((0:.05:1)');


parm.aggr.method='linear';
OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
% step = [0:.01:.1 .15:.05:.4 .5:.1:1];
% cst = (0:.05:1);
% [T1,T2]=meshgrid(step,cst);
% T0=zeros(numel(T1),4);
% T0(:,1)=T1(:);
% T0(:,2)=T1(:);
% T0(:,3)=T2(:);
% parm.n_real=48*size(T0,1);
%[d3.OF1, d3.OF2, d3.OF1std, d3.OF2std ] = OF(T0);



% OF1_range=[.019 .13];
% OF2_range=[3.2 7.6]*10^(-5);
out1 = (d.OF1-OF1_range(1))./(range(OF1_range));
out2 = (d.OF2-OF2_range(1))./(range(OF2_range));
out = out1+out2;
rout = d.OF1std./range(OF1_range) + d.OF2std./range(OF2_range);

figure(21);clf; hold on;
for i=1:numel(d.OF1)
    plot([d.OF1(i)-d.OF1std(i) d.OF1(i)+d.OF1std(i)],[d.OF2(i) d.OF2(i)],'k')
    plot([d.OF1(i) d.OF1(i)],[d.OF2(i)-d.OF2std(i) d.OF2(i)+d.OF2std(i)],'k')
end
% plot(d.OF1,d.OF2,'-k')
scatter(d.OF1,d.OF2,[],out,'filled')



figure(5); clf;hold on;
[~,min_out] = min(out);
subplot(1,3,1);hold on;
imagesc(T1(1,:),T2(:,1),reshape(d.OF1,numel(T1(1,:)),numel(T2(:,1))));
plot([T1(min_out) T1(min_out) 0],[0 T2(min_out) T2(min_out)],'r','linewidth',2);
xlabel('Step'); ylabel('cst');title('vario'); axis equal tight
subplot(1,3,2); hold on;
imagesc(T1(1,:),T2(:,1),reshape(d.OF2,numel(T1(1,:)),numel(T2(:,1)))); 
plot([T1(min_out) T1(min_out) 0],[0 T2(min_out) T2(min_out)],'r','linewidth',2);
xlabel('Step'); ylabel('cst');title('joint pdf'); axis equal tight
subplot(1,3,3); hold on;
imagesc(T1(1,:),T2(:,1),reshape(out,numel(T1(1,:)),numel(T2(:,1))));
plot([T1(min_out) T1(min_out) 0],[0 T2(min_out) T2(min_out)],'r','linewidth',2);
xlabel('Step'); ylabel('cst');title('combined'); axis equal tight


subplot(2,3,[4 6]); hold on;
my_col = repmat((out-min(out))/range(out),1,3); 
for i_t=1:numel(T1)
    plot([0 T1(i_t) T1(i_t) 1],[T2(i_t) T2(i_t) 0 0],'color',my_col(i_t,:));
end
plot([0 T1(min_out) T1(min_out) 1],[T2(min_out) T2(min_out) 0 0],'r','linewidth',2);
xlabel('Order of simulation'); ylabel('Weight value'); ylim([0 1])




figure(5); clf;
list = {d d01 d68};
for i=1:numel(list)
    d=list{i};
    out1 = (d.OF1-OF1_range(1))./(range(OF1_range));
    out2 = (d.OF2-OF2_range(1))./(range(OF2_range));
    out = out1+out2;
    subplot(1,3,1);hold on;
    imagesc(d.T1(1,:),d.T2(:,1),reshape(d.OF1,numel(d.T1(1,:)),numel(d.T2(:,1))));
    subplot(1,3,2); hold on;
    imagesc(d.T1(1,:),d.T2(:,1),reshape(d.OF2,numel(d.T1(1,:)),numel(d.T2(:,1))));
    subplot(1,3,3); hold on;
    imagesc(d.T1(1,:),d.T2(:,1),reshape(out,numel(d.T1(1,:)),numel(d.T2(:,1))));
end

[~,min_out] = min(out);
subplot(1,3,1);hold on;
plot([d.T1(min_out) d.T1(min_out) 0],[0 d.T2(min_out) d.T2(min_out)],'r','linewidth',2);
xlabel('Step'); ylabel('cst');title('vario'); axis equal tight
subplot(1,3,2); hold on;
plot([d.T1(min_out) d.T1(min_out) 0],[0 d.T2(min_out) d.T2(min_out)],'r','linewidth',2);
xlabel('Step'); ylabel('cst');title('joint pdf'); axis equal tight
subplot(1,3,3); hold on;
plot([d.T1(min_out) d.T1(min_out) 0],[0 d.T2(min_out) d.T2(min_out)],'r','linewidth',2);
xlabel('Step'); ylabel('cst');title('combined'); axis equal tight




parm.n_real=48*6;
OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
[~,~,Res,E1_x,E1_y,E2]=OF([.8 .4]);
figure(55);
subplot(1,2,1); imagesc(kern.axis_sec, kern.axis_prim, reshape(mean(E2,2),numel(kern.axis_prim),numel(kern.axis_sec))+kern.dens); 
xlim([2 12]); ylim([-10 -4]); caxis([0 max(kern.dens(:))]); set(gca,'Ydir','normal'); 
subplot(1,2,2); imagesc(kern.axis_sec, kern.axis_prim, kern.dens)
xlabel('\sigma [mS/m]'); ylabel('[Log m/s]'); xlim([2 12]); ylim([-10 -4]); caxis([0 max(kern.dens(:))]) ; set(gca,'Ydir','normal'); 

figure(56); hold on
addpath('/home/rnussbaumer/MATLAB/raacampbell_shadedErrorBar/raacampbell-shadedErrorBar-4c6d623/shadedErrorBar.m')
col=get(groot,'DefaultAxesColorOrder');

subplot(2,1,1); hold on;
shadedErrorBar(grid_gen.x(id.x),mean(E1_x(id.x,:),2)'+Gamma_id.x,[max(E1_x(id.x,:),[],2)-mean(E1_x(id.x,:),2)  abs(min(E1_x(id.x,:),[],2)-mean(E1_x(id.x,:),2))],{'color',col(1,:)},1);
h2(1)=plot(grid_gen.x(id.x),gamma_x_s(id.x),'k','linewidth',2);
h2(2)=plot(grid_gen.x(id.x),Gamma_id.x,'--k','linewidth',3);
axis tight;
subplot(2,1,2); hold on,
shadedErrorBar(grid_gen.y(id.y),mean(E1_y(id.y,:),2)'+Gamma_id.y,[max(E1_y(id.y,:),[],2)-mean(E1_y(id.y,:),2)  abs(min(E1_y(id.y,:),[],2)-mean(E1_y(id.y,:),2))],{'color',col(1,:)},1);
h2(1)=plot(grid_gen.y(id.y),gamma_y_s(id.y),'k','linewidth',2);
h2(2)=plot(grid_gen.y(id.y),Gamma_id.y,'--k','linewidth',3);
axis tight;


figure(57); hold on
caxis_v = [min(KlogNscore(:)) max(KlogNscore(:))];
subplot(4,1,1); imagesc(grid_gen.x,grid_gen.y,KlogNscore);
caxis(caxis_v)
subplot(4,1,2); imagesc(grid_gen.x,grid_gen.y, mean(Res,3))
caxis(caxis_v)
subplot(4,1,3); imagesc(grid_gen.x,grid_gen.y, Res(:,:,20))
caxis(caxis_v)
subplot(4,1,4); imagesc(grid_gen.x,grid_gen.y, std(Res,[],3))


figure(58); hold on
caxis_v = [min(log(K_true(:))) max(log(K_true(:)))];
subplot(2,1,1); imagesc(grid_gen.x,grid_gen.y, mean(Res,3))
%caxis(caxis_v)
subplot(2,1,2); imagesc(grid_gen.x,grid_gen.y,log(K_true));
caxis(caxis_v)


step = 0:.25:1;
T0=nan(numel(step),numel(step),numel(step),numel(step),4);
d4.OF1=nan(numel(step),numel(step),numel(step),numel(step));
d4.OF2=nan(numel(step),numel(step),numel(step),numel(step));
d4.OF1std=nan(numel(step),numel(step),numel(step),numel(step));
d4.OF2std=nan(numel(step),numel(step),numel(step),numel(step));

for i1=1:numel(step)
    for i2=1:numel(step)
        for i3=1:numel(step)
            for i4=1:numel(step)
                T0(i1,i2,i3,i4,:)=[step(i1) step(i2) step(i3) step(i4)];
                T0t=reshape(T0(i1,i2,i3,i4,:),1,4);
                [d4.OF1(i1,i2,i3,i4), d4.OF2(i1,i2,i3,i4), d4.OF1std(i1,i2,i3,i4), d4.OF2std(i1,i2,i3,i4) ] = OF(T0t);
            end
        end
    end
end

g=reshape(T0,625,4)';
OF1_range=[.019 .13];
OF2_range=[3.2 7.6]*10^(-5);
out = (reshape(d4.OF1,625,1)-OF1_range(1))./(range(OF1_range))+ (reshape(d4.OF2,625,1)-OF2_range(1))./(range(OF2_range));

figure(1002);
subplot(3,3,[1 3]);  histogram(out,200);
lim=[0 .3 .45 .7 1.3 1.9 Inf];
for i=1:numel(lim)-1
    subplot(3,3,3+i)
    id2 = out>lim(i) & out<lim(i+1);
    plot(g(:,id2),'k'); xlim([1 4]);ylim([0 1])
end


%% Classical Metropolis with Multi-grid weights. 4 params. which can be group. 

parm.aggr.method='hybrid';
parm.n_real=48*2;

n_pool=48;
parpool(n_pool);

nsamples = 60*50*4;
var_sampling = .05;
var_l = 1/50;

parm.aggr.step=[.01 .02 .03 .05 .08 .1 .2 .5];

T=nan(nsamples,numel(parm.aggr.step)+2);
out=nan(nsamples,1);
rout=nan(nsamples,1);
Acc = zeros(nsamples,1);

OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);


best=1;
T(1,:) = [1 1 1 1 1 1 0 0 0 1];

[~, out(1)]  = OF(T(1,:));

%i=1114
for i = 1165:nsamples

     T(i,:) = normrnd(T(best,:),var_sampling);
     T(i,1:end-1) = mod( T(i,1:end-1),T(i,end));

    disp(['Test nb ' num2str(i) ' with T=[' num2str(T(i,:)) ']'])
    [~, out(i)] = OF(T(i,:));
    

    % l = min( exp(-normrnd(L_m(i),L_s(i))/var_l)/exp(-normrnd(L_m(best), L_s(best))/var_l) ,1);
    l = min( exp(-out(i)/var_l)/exp(-out(best)/var_l) ,1);
    if rand <= l
        disp(['Accepted with out of : ' num2str(out(i)) ' | Best out ' num2str(out(best)) ' | likelihood ratio of ' num2str(l)])
        best=i;
        Acc(i) = 1;
    else
        disp(['Reject with out of : ' num2str(out(i)) ' Best out ' num2str(out(best)) ' | likelihood ratio of ' num2str(l)])
    end
%     
%     if mod(i,20)==0
%         if mean(Acc(i-19:2:i))>.5
%             var_sampling(2) = var_sampling(2)*(1+rand());
%             %var_l = var_l/(1+rand());
%         elseif mean(Acc(i-19:2:i))<.25
%             var_sampling(2) = var_sampling(2)/(1+rand());
%             %var_l = var_l*(1+rand());
%         end
%         if mean(Acc(i-18:2:i))>.5
%             var_sampling(1) = var_sampling(1)*(1+rand());
%             %var_l = var_l/(1+rand());
%         elseif mean(Acc(i-18:2:i))<.25
%             var_sampling(1) = var_sampling(1)/(1+rand());
%             %var_l = var_l*(1+rand());
%         end
%         disp(['Update var to: ' num2str(var_sampling) ' accpeted ratio of: ' num2str(mean(Acc(i-9:i)))])
%     end
end


figure(3); clf;
for i_parm=1:size(T,2)
    subplot(size(T,2),1,i_parm); hold on; histogram(T(Acc==1,i_parm)); 
end

figure(31); clf;
for i_parm=1:size(T,2)
    subplot(size(T,2),1,i_parm); hold on; histogram(T(:,i_parm)); 
end

figure(4); clf;
for i_parm=1:size(T,2)
    subplot(size(T,2),1,i_parm); hold on; scatter(T(:,i_parm),out); ylim([0 .30]);xlim([0 1.6])
end

figure(5); hold on;
cmap = colormap; c = round(1+(size(cmap,1)-1)*(out - min(out(Acc>0)))/(max(out(Acc>0))-min(out(Acc>0))));
for i_sim=1:i
    if Acc(i_sim)>0
     plot3([0 parm.aggr.step]+diff([0 parm.aggr.step 1])/2, T(i_sim,1:end-1)*T(i_sim,end),repmat(-out(i_sim),numel(parm.aggr.step)+1,1),'color',cmap(c(i_sim),:))
    end
end

figure(6); hold on
[~,idx]=min(out);
stairs([0 parm.aggr.step 1], [T(idx,1:end-1) T(idx,end-1)]*T(idx,end))
stairs([0 parm.aggr.step 1], (1-[T(idx,1:end-1) T(idx,end-1)])*T(idx,end))


[~,idx]=min(out);
RestB = OF(T(idx,:));
figure(7);
for i_t=1:5
    subplot(6,1,i_t); 
    r= RestB(:,:,i_t);
    imagesc(grid_gen.x,grid_gen.y,inve((r-mean(r(:)))/std(r(:))));
    caxis([kern.axis_prim(1) kern.axis_prim(end)])
end
subplot(6,1,1); 
imagesc(grid_gen.x,grid_gen.y,std(RestB,[],3));


figure(5); clf;hold on;
F= scatteredInterpolant(T(1:i,1),T(1:i,2),L_m(1:i));
[Xq,Yq] = meshgrid(0:.01:1,0:.01:1);
imagesc(Xq(1,:),Yq(:,1)',F(Xq,Yq));
scatter(T(1:i,1),T(1:i,2),L_m(1:i),'o','filled')

figure(5); hold on;
my_col = repmat((L_m(Acc==1)-min(L_m))/range(L_m(Acc==1)),1,3); 
Tt=T(Acc==1,:);
for i_t=1:size(Tt,1)
    plot([0 Tt(i_t,1) Tt(i_t,1) 1],[Tt(i_t,4) Tt(i_t,4) Tt(i_t,5) Tt(i_t,5)],'color',my_col(i_t,:));
end




%% Personal Optimization procedure, based on MH but a bit different

        
        
parm.n_real=4;
var_sampling=1;

%n_pool=48;
%parpool(n_pool);

n_parm=10;

T=cell(n_parm,1);
OF=cell(n_parm,1);
OFr=cell(n_parm,1);
A=cell(n_parm,1);

n_last_iter=6;


for i1=1:n_parm % number of param to calibrate
    
    T{i1} = cell(i1,1);
    OF{i1} = cell(i1,1);
    OFr{i1} = cell(i1,1);
    A{i1} = cell(i1,1);
    
    if i1==1
        T{i1}{1} = .5 ;
        [OF{i1}{1}, OFr{i1}{1}] = fmin(T{i1}{1},parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
        A{i1}{1} = 1 ;
        disp('Initialized the first test case')
    else
        T{i1}{1} = T{i1-1}{end}(end); % upodate the size of T_best
        disp(['Changed to new level with ' num2str(i1) ' params'])
    end
    
    for i2 = 1:i1 % each of param to calibrate
        
        if i2~=1
            disp(['Changed to the ' num2str(i2) 'th parameters of the ' num2str(i1) ' from this level'])
            T{i1}{i2} = T{i1}{i2-1}(end);
            OF{i1}{i2} = OF{i1}{i2-1}(end);
            OFr{i1}{i2} = OFr{i1}{i2-1}(end);
            A{i1}{i2} = 1;
        end

        condition=1;
        
        last=1;
        cur=2;
        i_var=0;
        var_sampling = exp(i_var);
        while condition 
            disp(['New test nb ' num2str(cur) ' | parameters: ' num2str(i2) '/' num2str(i1)])
            if (mod(cur,n_last_iter)==0)
                dacc_ratio = mean(A{i1}{i2}(cur-n_last_iter+1:cur-1));
                if dacc_ratio > .5 % accept too much -> reduce variance
                    i_var = i_var-1; 
                elseif dacc_ratio < .1 % reject too much -> increase variance
                    i_var = i_var-1; %
                end
                var_sampling = exp(i_var);
                disp(['dacc_ratio=' num2str(dacc_ratio) ' | var=' num2str(var_sampling)])
                if i_var<-5
                    condition=0;
                end
            end
            
            % Perturbation of the ith param with a var var and block
            % between 0 and 1
            T{i1}{i2}(cur,:) = mod(normrnd(T{i1}{i2}(last,:),var_sampling),1);
            disp(['New T=' num2str(T{i1}{i2}(cur,:)) ' | Last T =' num2str(T{i1}{i2}(last,:))])
            

            [OF{i1}{i2}(cur), OFr{i1}{i2}(cur)] = fmin(T{i1}{i2}(cur,:),parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
            
            while abs(OF{i1}{i2}(cur) - OF{i1}{i2}(last)) < abs(OFr{i1}{i2}(cur))+abs(OFr{i1}{i2}(last))
                 parm.n_real = parm.n_real*2;
                 disp(['Update the number of real to: ' num2str(parm.n_real)])
                 if parm.n_real >500
                    condition = 0;
                    disp('Too many real, stop the current parameter');
                    continue
                 end
                 [OF{i1}{i2}(last),OFr{i1}{i2}(last)] = fmin(T{i1}{i2}(last,:),parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
                 [OF{i1}{i2}(cur),OFr{i1}{i2}(cur)] = fmin(T{i1}{i2}(cur,:),parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd);
            end
            
            
            if OF{i1}{i2}(cur) < OF{i1}{i2}(last)
                last = cur;
                A{i1}{i2}(cur)=1;
                disp(['Accepted with OF=' num2str(OF{i1}{i2}(cur)) ' +/- ' num2str(OFr{i1}{i2}(cur))]);
            else
                A{i1}{i2}(cur)=0;
            end
            
            
%             if sum(A{i1}{i2}) >= n_last_iter
%                 y=OF{i1}{i2}(A{i1}{i2});
%                 y=y(en-n_last_iter:end);
%                 x=(1:n_last_iter)';
%                 slope = x'\y;
%                 ttest = ( (slope-0)*sqrt(n_last_iter-2) )/sqrt( sum((y-slope).^2) / sum( ( x-mean(x) ).^2 ));
%                 
%                 if (ttest< threasold)
%                     condition=0;
%                 end
%             end
            
            cur=cur+1;
        end
        
    end

end




















