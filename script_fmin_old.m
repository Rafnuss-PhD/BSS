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
parm.k.covar = kriginginitiaite(parm.k.covar);

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
id = grid_gen.x<parm.k.covar(1).range0(2).*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
Gamma_t_id = Gamma_t(id);
XY = kern.XY;
XY(:,2)= Nscore.forward(XY(:,2));
Sigma_d = Sigma.d(:);
dens = kern.dens(:)./sum(kern.dens(:));


%% Run Fminc
parm.n_real=48*6;

parm.aggr.method='linear';

OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);

OF([2]');


T0 = [0.0097    0.1019   0.7863 0.3586 ]; % T_1, T_2 w_1 w_2,
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
'PlotFcn',{@optimplotx, @optimplotfval , @optimplotstepsize  },'UseParallel',true);
x = fmincon(OF,T0,A,b,[],[],lb,ub,[],options)


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

%% Metroplis H

parm.aggr.method='mg';
parm.n_real=48*6;
OF = @(T) fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);

n_pool=48;
parpool(n_pool);

nsamples = 14*60;
T=cell(n_pool,1);
L=cell(n_pool,1);
Acc=cell(n_pool,1);
U=cell(n_pool,1);

parfor i_pool=1:n_pool

    T{i_pool}=nan(nsamples,sn);
    T{i_pool}(1,:) = unifrnd(0,1,1,sn);
    L{i_pool}=nan(nsamples,1);
    Acc{i_pool} = zeros(nsamples,1);
    U{i_pool}=rand(nsamples);
    
    best=1;
    T_b=T{i_pool}(1,:);
    L_b=OF(T_b);

    for i = 1:nsamples
        T{i_pool}(i,:) = max(0,min(1,normrnd(T_b,.1)));
        L{i_pool}(i) = OF(T{i_pool}(i,:));
        if U{i_pool}(i)<=min(L{i_pool}(i)/L_b,1)
            best=i;
            Acc{i_pool}(i) = 1;
        end
    end
    
    
end

t=reshape([T{:}],size(T{1},1),size(T{1},2),n_pool);
l=[L{:}];
acc=[Acc{:}];


save('result-BSS/MH_1.mat','nsamples','T','L','parm','Acc')

figure(1);clf; hold on;
for i_pool=1:n_pool
    plot(mean(T{i_pool}(Acc{i_pool}==1,:)))
end

for i_pool=1:min(5,n_pool)
    %subplot(2,3,i_pool); hold on;
    %plot(T{i_pool}(Acc{i_pool}==1,:)')
    %plot(mean(T{i_pool}(Acc{i_pool}==1,:)),'--k','linewidth',2)
    %subplot(2,3,6); hold on;
    plot(mean(T{i_pool}(Acc{i_pool}==1,:)))
end

figure(2); clf;
thr=.4;
for i_s=1:sn
    subplot(2,5,i_s); hold on; xlabel(num2str(i_s))
    tis = reshape(t(:,i_s,:),nsamples*n_pool,1);
    scatter(tis(l(:)>thr),l(l(:)>thr),'.k')
end


figure(3); clf;
for i_s=1:sn
    subplot(2,5,i_s); hold on; xlabel(num2str(i_s))
    tis = reshape(t(:,i_s,:),nsamples*n_pool,1);
    histogram(tis(acc(:)==1))
    axis tight;
end

%% Metroplis 2
        
        
parm.n_real=4;
var=1;

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
        [OF{i1}{1}, OFr{i1}{1}] = fmin(T{i1}{1},parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);
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
        var = exp(i_var);
        while condition 
            disp(['New test nb ' num2str(cur) ' | parameters: ' num2str(i2) '/' num2str(i1)])
            if (mod(cur,n_last_iter)==0)
                dacc_ratio = mean(A{i1}{i2}(cur-n_last_iter+1:cur-1));
                if dacc_ratio > .5 % accept too much -> reduce variance
                    i_var = i_var-1; 
                elseif dacc_ratio < .1 % reject too much -> increase variance
                    i_var = i_var-1; %
                end
                var = exp(i_var);
                disp(['dacc_ratio=' num2str(dacc_ratio) ' | var=' num2str(var)])
                if i_var<-5
                    condition=0;
                end
            end
            
            % Perturbation of the ith param with a var var and block
            % between 0 and 1
            T{i1}{i2}(cur,:) = mod(normrnd(T{i1}{i2}(last,:),var),1);
            disp(['New T=' num2str(T{i1}{i2}(cur,:)) ' | Last T =' num2str(T{i1}{i2}(last,:))])
            

            [OF{i1}{i2}(cur), OFr{i1}{i2}(cur)] = fmin(T{i1}{i2}(cur,:),parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);
            
            while abs(OF{i1}{i2}(cur) - OF{i1}{i2}(last)) < abs(OFr{i1}{i2}(cur))+abs(OFr{i1}{i2}(last))
                 parm.n_real = parm.n_real*2;
                 disp(['Update the number of real to: ' num2str(parm.n_real)])
                 if parm.n_real >500
                    condition = 0;
                    disp('Too many real, stop the current parameter');
                    continue
                 end
                 [OF{i1}{i2}(last),OFr{i1}{i2}(last)] = fmin(T{i1}{i2}(last,:),parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);
                 [OF{i1}{i2}(cur),OFr{i1}{i2}(cur)] = fmin(T{i1}{i2}(cur,:),parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd);
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




















