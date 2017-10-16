function [Rest,ww, t, k, parm] = BSS(nx,ny,hd,f0,sec,parm)

tik.global = tic;
%% * *INPUT CEHCKING*
% This section of the code generates a valid parm structure based on the
% inputted parm. If some parm value haven't been input, this section will
% fill automatically with defautl value. This may not allowed be the best.


% Paramter settings
if ~isfield(parm, 'seed_path'),     parm.seed_path      = 'shuffle'; end
if ~isfield(parm, 'seed_U'),        parm.seed_U         = 'shuffle'; end
if ~isfield(parm, 'seed_search'),   parm.seed_search    = 'shuffle'; end
if ~isfield(parm, 'saveit'),        parm.saveit         = 0; end % bolean, save or not the result of simulation
if ~isfield(parm, 'name'),          parm.name           = ''; end % name use for saving file
if ~isfield(parm, 'n_real'),        parm.n_real         = 1; end

% Kriging parameter
parm.k.covar = kriginginitiaite(parm.k.covar);
if ~isfield(parm, 'k') || ~isfield(parm.k, 'method'),   parm.k.method = 'sbss'; end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'lookup'),   parm.k.lookup = false; end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'nb'),       parm.k.nb = 30; end

if ~isfield(parm, 'k') || ~isfield(parm.k, 'wradius')
    parm.k.wradius  = 3;
end
k = parm.k;

% Path
if ~isfield(parm, 'path'),          parm.path            = 'linear'; end
if ~isfield(parm, 'path_random'),   parm.path_random     = 1; end
if ~isfield(parm, 'mg'),            parm.mg              = 1; end

% Hard data
hd.x=hd.x(:); hd.y=hd.y(:); hd.d=hd.d(:); 
if ~isfield(hd, 'n'),            hd.n              = numel(hd.d); end
if ~isfield(hd, 'id'),          
    hd.id = sub2ind([ny nx],hd.y,hd.x);
else
    hd.id=hd.id(:);
end
hd.scale=nan(hd.n,1);

% Secondary variable
assert(all(size(sec.pdf)==[numel(sec.axis) nx*ny]), 'sec.pdf in not ok')

% assert(size(parm.aggr.T,2)<3)

%% 1. Creation of the grid an path
[Y, X] = ndgrid(1:ny,1:nx);


%% 2. Define Path
tik.path = tic;
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
t.path = toc(tik.path);


%% 3. Initialization Spiral Search

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

%% 3. Initialization Covariance Lookup Table
if k.lookup
    ss_a0_C = zeros(ss_n,1);
    ss_ab_C = zeros(ss_n);
    for i=1:numel(k.covar)
        a0_h = sqrt(sum(([ss_Y_s(:) ss_X_s(:)]*k.covar(i).cx).^2,2));
        ab_h = squareform(pdist([ss_Y_s ss_X_s]*k.covar(i).cx));
        
        ss_a0_C = ss_a0_C + kron(k.covar(i).g(a0_h), k.covar(i).c0);
        ss_ab_C = ss_ab_C + kron(k.covar(i).g(ab_h), k.covar(i).c0);
    end
end
% Transform ss.ab_C sparse?



%% 3. Simulation
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
    if k.lookup
        ss_a0_C_s = ss_a0_C(ss_id);
        ss_ab_C_s = ss_ab_C(ss_id,ss_id);
    end
    
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
            
            if k.lookup
                a0_C = ss_a0_C_s(neigh(1:n));
                ab_C = ss_ab_C_s(neigh(1:n), neigh(1:n));
            else
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

if parm.saveit
    filename=['result-SGS/SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    mkdir('result-SGS/')
    save(filename, 'parm', 'nx','ny','start','nb', 'path', 'sn', 'k','NEIGH','S','LAMBDA')
end

%% Realization loop

tik.real = tic;
Rest = nan(ny,nx,parm.n_real);
parm_seed_U = parm.seed_U;


for i_real=1:parm.n_real
    Res=nan(ny,nx);
    Res(hd.id) = hd.d;
    rng(parm_seed_U);
    U=rand(ny,nx);
    for i_scale = 1:sn
        for i_pt = start(i_scale)+(1:nb(i_scale))
            n = ~isnan(NEIGH(i_pt,:));
            fkrig = normpdf(sec.axis, LAMBDA(i_pt,n)*Res(NEIGH(i_pt,n))', sqrt(S(i_pt)));
            fkrig = fkrig' ./ sum(fkrig);
            
            fsec = sec.pdf(:,path(i_pt));
            
            w = aggr_fx(i_real,i_pt/sum(nb),parm.aggr,i_scale/sn);
            
            
            fa = f0.^0 .* fkrig.^(1-w) .* fsec.^w;
            
            cfa = cumsum([0 ; fa(2:end-1)+eps ; eps]) ./ (sum(fa(2:end-1)+eps));
            Res(path(i_pt)) =  interp1(cfa, sec.axis, U(i_pt),'linear');
            
           
            if 0==1
                figure(3); clf; hold on;
                plot(sec.axis,f0)
                plot(sec.axis,fsec)
                plot(sec.axis,fkrig)
                plot(sec.axis,fa./sum(fa))
                legend('prior', 'sec', 'krig', 'aggr')
                
                figure(1); clf;
                imagesc(Res); hold on;
                plot(X(path(i_pt)), Y(path(i_pt)),'or')
                plot(X(NEIGH(i_pt,n)), Y(NEIGH(i_pt,n)),'xk')
                caxis([-4 4]); axis equal
                keyboard
            end
        end
    end
    Rest(:,:,i_real) = Res;
end



if parm.saveit
    filename=['result-SGS/SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    mkdir('result-SGS/')
    save(filename, 'parm','nx','ny', 'Rest', 't', 'k','U')
end

t.real = toc(tik.real);
t.global = toc(tik.global);
disp(['Run finished in ' num2str(t.global*60)] )
end


function w = aggr_fx(i_real, x, aggr, s)


i_t = mod(i_real,size(aggr.T,1));
if i_t==0;
    i_t=size(aggr.T,1);
end
assert(x<=1,'error')

switch aggr.method
    case 'cst'
        w=aggr.T(i_t);
    case 'step'
        if (x<aggr.T(i_t))
            w=0;
        else 
            w=1;
        end
    case 'linear'
        if (x<aggr.T(i_t,1))
            w  = aggr.T(i_t,3);
        elseif (x>aggr.T(i_t,2))
            w = aggr.T(i_t,4);
        else 
            w =  aggr.T(i_t,3) + ( x - aggr.T(i_t,1) )/(aggr.T(i_t,2)-aggr.T(i_t,1)) * (aggr.T(i_t,4)-aggr.T(i_t,3));
        end
    case 'sigmoid'
        a = aggr.T(i_t,1);
        b = aggr.T(i_t,2);
        w = (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
    case 'mg'
        i_s = ceil(s*size(aggr.T,2));
        w = aggr.T(i_t,i_s);
    otherwise
        error('no Aggr method defined')  
end
end
