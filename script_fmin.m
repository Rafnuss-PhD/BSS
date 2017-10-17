%% Preparation
load('result-BSS/GEN-Run_1_2017-05-07_14-37');
addpath(genpath('./.'));
Nscore = nscore(kern, struct('nscore', 1), 0); %Prim.d, kern.axis_prim, 'pchip', 'pchip', parm.plot.ns
[~,s_id]=min(bsxfun(@minus,kern.axis_sec,Sigma.d(:)).^2,[],2);
sec_pdf = kern.dens(:,s_id);
sec.pdf = bsxfun(@times, sec_pdf, 1./sum(sec_pdf));
sec.axis = Nscore.forward(kern.axis_prim);

parm.k.covar = gen.covar;
parm.k.covar.range0 = fliplr(gen.covar.range0) ./ [grid_gen.dy grid_gen.dx];

parm.k.wradius = 1.3;
parm.k.lookup = false;
parm.k.nb = 40;
parm.path_random     = 1;
parm.path            = 'linear';
parm.mg=1;

parm.seed_path      = 'shuffle';
parm.seed_search    = 'shuffle';

parm.aggr.method='linear';
parm.n_real=1;

% use the log of hyd. cond.
hd = sampling_pt(struct('x',1:grid_gen.nx,'y',1:grid_gen.ny),K_true,2,0);
hd.d = Nscore.forward(hd.d);
hd.scale=nan(hd.n,1);
hd.x=hd.x(:); hd.y=hd.y(:); hd.d=hd.d(:); 

f0=kern.prior ./ sum(kern.prior);
nx = grid_gen.nx;
ny = grid_gen.ny;

k = parm.k;

%% BSS function

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
x = ceil( min(k.covar(1).range(2)*k.wradius, nx));
y = ceil( min(k.covar(1).range(1)*k.wradius, ny));
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


%% Prepare OF
id = grid_gen.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
Gamma_t_id = Gamma_t(id);
XY = kern.XY;
XY(:,2)= Nscore.forward(XY(:,2));
Sigma_d = Sigma.d(:);
dens = kern.dens(:)./sum(kern.dens(:));


%% Run Fminc
parm.n_real=48;

parm.aggr.method='linear';

OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);


T0 = [0.0097    0.1019   0.7863 0.3586 ]; % T_1, T_2 w_1 w_2,
T0 = [0.99    0.991   0.999 0.62 ]; % T_1, T_2 w_1 w_2,
T0 = [0    0    .5    .5]; % T_1, T_2 w_min w_max,
out = OF(T0);

out = OF(T0);
out = OF([3 NaN 0 NaN]);

A=[1 -1 0 0 ; 0 0 -1 1];
b=[0;0];
lb=[0 0 0 0];
ub=[1 1 1 1];

options = optimoptions('fmincon','Display','iter-detailed','Diagnostics','on','MaxFunctionEvaluations',100,...
'PlotFcn',{@optimplotx, @optimplotfval , @optimplotstepsize  });
x = fmincon(OF,T0,A,b,[],[],lb,ub,[],options)


options = optimoptions(@particleswarm,'PlotFcn',,'Display','iter','UseVectorized',true);
x = particleswarm(OF,4,lb,ub);


options = optimoptions(@simulannealbnd,'PlotFcn',{@saplotbestf, @saplotbestx, @saplotf, @saplotx, @saplotstopping , @saplottemperature},'Display','iter','PlotInterval',10);
x = simulannealbnd(OF,T0,lb ,ub,options);


options = optimoptions(@patternsearch,'PlotFcn',{@psplotbestf, @psplotmeshsize, @psplotfuncount, @psplotbestx},'Display','iter','PlotInterval',1);
x = patternsearch(OF,T0,A,b,[],[],lb,ub,[],options);
