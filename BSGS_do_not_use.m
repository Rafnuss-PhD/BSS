%% BSGS is computing the Bayesian Sequential Gaussian Simulation.
% BSGS is a sequential stochatic simulation of a primary variable from
% combination of a primary variable (X) and a secondary data (Z). 
% The primary variable can be at any point of this grid (X.x,X.y)
%
% INPUT:
%
% * X.d     : Primary variable (1D-X.n elts)
% * X.x,y   : x,y-coordinate of the primary variable (1D-X.n elts)
% * Z.d     : Secondary variable (2D- grid.nx x grid.n elts)
% * Z.x,y   : x,y-coordinate of the secondary variable (1D elts)
% * Z.std   : Secondary variable std error (2D-grid.nx x grid.n elts)
% * grid    : grid informations
% * parm    : Parameters of the simulation.
%   * gen ([]): generation parmater structure (ouput of data_generation.m)
%   * likelihood (1): bolean. with or without using likelihood.
%   * scale (1:grid{end}): array of the scales (i, where nx=2^grid{i}+1) to simulate 
%   * name (''): name of the simulation
%   * seed(rand): random number. use to reproduce exactly the same simulation.
%   * saveit (1): save or not in a file
%   * unit (''): unit of the primary variable, used for plot
%   * n_realisation (1): number of realisations.
%   * neigh (1): bolean. with or without "smart neighbour"
%   * nscore (1): bolean. with or without  normal score transform
%   * cstk (1): bolean. with or without constant path
%   * cstk_s (0 if cstk, Inf else): scale at which cst path is switch on
%   * plot:
%       * bsgs (0): plot. ?
%       * ns (0): plot. ?
%       * sb (0): plot. ?
%       * kernel (0): plot. ?
%       * fitvar (0): plot. ?
%       * krig (0): plot. ?
%   * kernel_range: kernel create a grid of X and Z where the range of the grid
%       is define with min and max of X and Z
%       * min (ceil(grid{end}.x(end)/parm.k.model(1,2)*3)): 
%       * max (ceil(grid{end}.x(end)/parm.k.model(1,2)*3)):
%   * k:
%       * nb_neigh: row: min and max number of point to use in each
%       quadrant (first 4 columns) and hard data (5th col.)
%       * model: covariance model.
%       * var: variance of each variogram
%       * sb:
%           * nx:
%           * ny:
%       * wradius (1.3): 
% 
%
% OUTPUT:
%
% * Y       : Simulated Primary variable
% * t       : Time of simulations
% * kernel  : kernel information
% * k       : kriging information
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% Referances:
%

function [Y, t, kernel, k] = BSGS(X,Z,X_true,grid_gen,parm)
t.tic.global = tic;
addpath(genpath('./.'))

%% * *INPUT CEHCKING*
% force input in column
X.x=X.x(:);X.y=X.y(:);X.d=X.d(:);Z.x=Z.x(:);Z.y=Z.y(:);

% default value
if ~isfield(parm, 'seed'),          parm.seed           = rand(); end
if ~isfield(parm, 'saveit'),        parm.saveit         = 1; end % bolean, save or not the result of simulation
if ~isfield(parm, 'name'),          parm.name           = ''; end % name use for saving file
if ~isfield(parm, 'familyname'),    parm.familyname     = ''; end
if ~isfield(parm, 'unit'),          parm.unit           = ''; end % unit of the primary variable, used in plot 
if ~isfield(parm, 'n_realisation'), parm.n_realisation  = 1; end
if ~isfield(parm, 'scale')
    parm.scale = repmat(1:max([grid_gen.sx,grid_gen.sy]),2,1); 
    parm.scale(1,parm.scale(1,:)>grid_gen.sx) = grid_gen.sx;
    parm.scale(2,parm.scale(2,:)>grid_gen.sy) = grid_gen.sy;
end
if ~isfield(parm, 'p_w'),
    parm.p_w     = repmat(0.5, 1, size(parm.scale,2)); 
elseif numel(parm.p_w)==1
    parm.p_w     = repmat(parm.p_w, 1,size(parm.scale,2)); 
elseif numel(parm.p_w) == size(parm.scale,2)
    parm.p_w    = parm.p_w;
else
    parm.p_w     = linspace(parm.p_w(1), parm.p_w(end),size(parm.scale,2)); 
end
% Run option
if ~isfield(parm, 'neigh'),         parm.neigh          =1; end % smart-neighbouring activated or not
if ~isfield(parm, 'nscore'),        parm.nscore         =1; end % use normal score (strongly advice to use it.)
if ~isfield(parm, 'cstk_s')% cstk_s is the scale at which cst is switch on
    if ~isfield(parm, 'cstk'),      parm.cstk           = 1; end % constant path and kriging weight activated or not
    if parm.cstk
        parm.cstk_s = 0; % will never use cstk
    else
        parm.cstk_s = Inf; % will always use cstk
%     end
end
if ~isfield(parm, 'fitvar'), parm.fitvar     =0; end % fit the variogram to the data or used the given one in parm.covar
% Plot
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'bsgs'),  parm.plot.bsgs   =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'ns'),    parm.plot.ns     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'sb'),    parm.plot.sb     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'kernel'),parm.plot.kernel =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'fitvar'),parm.plot.fitvar =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'krig'),  parm.plot.krig   =0; end
% Kernel parameter
if ~isfield(parm, 'kernel_range') || ~isfield(parm.kernel_range, 'min') || ~isfield(parm.kernel_range, 'max')% range used in the kernel estimator
    parm.kernel_range.min = [min(Z.d(:))-.2*range(Z.d(:)) min(X.d(:))-.2*range(X.d(:))];
    parm.kernel_range.max = [max(Z.d(:))+.2*range(Z.d(:)) max(X.d(:))+.2*range(X.d(:))];
end
% Kriging parameter
if ~isfield(parm, 'k') || ~isfield(parm.k, 'nb_neigh')
    parm.k.nb_neigh = [0 0 0 0 0; 5 5 5 5 5];
end 
if ~isfield(parm, 'k') || ~isfield(parm.k, 'model') || ~isfield(parm.k, 'var')
    assert(isfield(parm, 'gen'),'You need to define parm.covar or parm.gen')
    parm.k.model    = parm.gen.covar.modele; 
    parm.k.var      = parm.gen.covar.c; 
end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'sb') || ~isfield(parm.k.sb, 'nx') || ~isfield(parm.k.sb, 'ny') % super-block grid size (hard data)
    parm.k.sb.nx    = ceil(grid_gen.x(end)/parm.k.model(1,2)*3);
    parm.k.sb.ny    = ceil(grid_gen.y(end)/parm.k.model(1,3)*3);
end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'wradius')
    parm.k.wradius  = 1.3; 
end
parm.k.rotation     = parm.k.model(1,4);
parm.k.range        = parm.k.model(1,2:3);


if parm.fitvar
    k = fit_variogramm(X,Z,parm, X_true);
else
    k = parm.k;
end

% Compute the rot matrix
for i=1:size(parm.k.model,1)
    ang=parm.k.model(i,4); cang=cos(ang/180*pi); sang=sin(ang/180*pi);
    rot = [cang,-sang;sang,cang];
    parm.k.cx{i} = rot/diag(parm.k.model(i,2:3));
end

% Check the input for correct size dans dimension
assert(ismatrix(Z.d),'Z is not 2Y.mat_simD');
assert(all([numel(Z.y), numel(Z.x)]==size(Z.d)),'Z.x,y does not seems to be ok');
assert(size(X.d,2)<=1,'X.d is not a vertical vector 1D');
assert(size(X.x,2)<=1,'X.dx is not a vertical vector 1D');
assert(size(X.y,2)<=1,'X.y is not a vertical vector 1D');
assert(all(size(X.y)==size(X.x)),'X.x and X.y don''t have the same dimension');
assert(all(size(X.d)==size(X.x)),'X.d and X.x (or X.y) don''t have the same dimension');
assert(max(parm.scale(:))<=max([grid_gen.sx,grid_gen.sy]),'This scale of simulation does not exist')
assert(all(parm.scale(1,:)<=grid_gen.sx),'monotonicly increasing scale')
assert(all(parm.scale(2,:)<=grid_gen.sy),'monotonicly increasing scale')

% creation of the grid
parm.n_scale=size(parm.scale,2);
grid=cell(parm.n_scale,1);
for scale_i = 1:parm.n_scale
    grid{scale_i}.sx=parm.scale(1,scale_i);
    grid{scale_i}.sy=parm.scale(2,scale_i);
    grid{scale_i}.nx=2^grid{scale_i}.sx+1;
    grid{scale_i}.ny=2^grid{scale_i}.sy+1;
    grid{scale_i}.nxy=grid{scale_i}.nx*grid{scale_i}.ny; % total number of cells

    grid{scale_i}.dx=grid_gen.x(end)/(grid{scale_i}.nx-1);
    grid{scale_i}.dy=grid_gen.y(end)/(grid{scale_i}.ny-1);

    grid{scale_i}.x=linspace(0, grid_gen.x(end), grid{scale_i}.nx); % coordinate of cells center
    grid{scale_i}.y=linspace(0, grid_gen.y(end), grid{scale_i}.ny);
    grid{scale_i}.xy=1:grid{scale_i}.nxy;

    [grid{scale_i}.X, grid{scale_i}.Y] = meshgrid(grid{scale_i}.x,grid{scale_i}.y); % matrix coordinate
end



%% * 2. *NON-PARAMETRIC RELATIONSHIP*
% The joint pdf of the primary and secondary is build using a bivariate 
% kernel density estimator (Botev, Grotowski, & Kroese, 2010).
kernel = kernel_est(X, Z, parm.kernel_range, parm.plot.kernel);


%% * 3. *NORMAL SCORE TRANSFORM*
% Based on the hard data (well samples), a normal score transform is 
% created using interpolation with power or exponential tail extrapolation.
% Matlab symbolique function are used for efficient coding. The back 
% transform of the prior normal distribution function is also created 
% (return the pdf in the initial space from the mean and variance in the 
% normal space)
if parm.nscore
    Nscore = nscore_perso(X.d, 'linear', 'linear', kernel, parm.plot.ns);
else
    Nscore.forward = @(x) x;
    Nscore.inverse = @(x) x;
    Nscore.dist    = @(mu,sigma) normpdf(kernel.y,mu,sigma)/sum(normpdf(kernel.y,mu,sigma)); 
end


%% * 1. *SUPERBLOCK GRID CREATION* 
% A mask (Boolean value) of the hard data is assigned to each superblock 
% as follow: Only the n-closest (normalized by the covariance range) points
% (inside the ellipse/windows) to the centre of the superblock will be 
% true. During the kriging, the mask of the superblock of the estimated 
% point will be used to select the hard to add to the system
if parm.neigh
    k.sb.nx = parm.k.sb.nx; % number of superblock grid
    k.sb.ny = parm.k.sb.ny;
    [k, X] = SuperBlockGridCreation(k, grid_gen.x(end), grid_gen.y(end), X, parm.plot.sb, parm.k.nb_neigh(2,5));
end



%% * 4. *SIMULATION*
% Create the normal space primary variable of known data
X.d_ns = Nscore.forward(X.d);
Y=cell(parm.n_scale,1); % allocate variable space
X_ini = X; %keep initial data intact
if parm.neigh; k.sb.mask_ini = k.sb.mask; end

for scale_i=1:parm.n_scale % for each scale
    t.tic.scale = tic;
    
    %%
    % * *INITIATE SCALE SIMULATION*
    % Allocating space for resulting field.
    Y{scale_i}.x=grid{scale_i}.x; 
    Y{scale_i}.y=grid{scale_i}.y; 
    Y{scale_i}.X=grid{scale_i}.X; 
    Y{scale_i}.Y=grid{scale_i}.Y; 
    Y{scale_i}.nx=grid{scale_i}.nx; 
    Y{scale_i}.ny=grid{scale_i}.ny;
    Y{scale_i}.m=repmat({nan(grid{scale_i}.ny,grid{scale_i}.nx)},parm.n_realisation,1); % matrix des resutlats
    
    % Populate the grid from previous scale.
    if scale_i~=1 % not at first scale
        for i_sim=1:parm.n_realisation
            Y{scale_i}.m{i_sim}( 1:(grid{scale_i-1}.dy/grid{scale_i}.dy):end, 1:(grid{scale_i-1}.dx/grid{scale_i}.dx):end) = Y{scale_i-1}.m{i_sim};
        end
    end
    
    % Assimilate the hard data (X) into the grid
    hard_data_idx=find(ismember(X.y,grid{scale_i}.Y)&ismember(X.x,grid{scale_i}.X));
    for i_sim=1:parm.n_realisation
        for hd=1:numel(hard_data_idx)
            Y{scale_i}.m{i_sim}(X.x(hard_data_idx(hd))==grid{scale_i}.X & X.y(hard_data_idx(hd))==grid{scale_i}.Y) = X.d(hard_data_idx(hd));
        end
    end
    % Remove the assimilated data from X.
    X.d(hard_data_idx)=[];
    X.x(hard_data_idx)=[];
    X.y(hard_data_idx)=[];
    X.d_ns(hard_data_idx)=[];
    if parm.neigh; k.sb.mask(:,:,hard_data_idx)=[]; end
    X.n=numel(X.d);
    
    % Create the normal space result matrix and populate with known value 
    Y{scale_i}.m_ns = repmat({nan(size(Y{scale_i}.m{1}))},parm.n_realisation,1); 
    I=~isnan(Y{scale_i}.m{1});
    for i_realisation=1:parm.n_realisation
        Y{scale_i}.m_ns{i_realisation}(I) = Nscore.forward(Y{scale_i}.m{i_realisation}(I));
    end


    
    %%
    % * *SPIRAL SEARCH*: create the windows for kringing with the function to compute
    % the normalized distence and the order of visit of the cells.
    % Spiral Search setting: previously data (on grid{scale_i} location)
    if parm.neigh
        [k.ss.el.X, k.ss.el.Y] = meshgrid(0:max(ceil(k.range(1)*k.wradius/grid{scale_i}.dx),ceil(k.range(2)*k.wradius/grid{scale_i}.dy)));% grid{scale_i} of searching windows
        [k.ss.el.X_T, k.ss.el.Y_T]=rotredtrans(k.ss.el.X*grid{scale_i}.dx, k.ss.el.Y*grid{scale_i}.dy, k.rotation, k.range); % transforms the grid{scale_i}
        k.ss.el.dist = sqrt(k.ss.el.X_T.^2 + k.ss.el.Y_T.^2); % find distence
        [k.ss.el.dist_s, k.ss.el.dist_idx] = sort(k.ss.el.dist(:)); % sort according distence.
        k.ss.el.X_s=k.ss.el.X(k.ss.el.dist_idx); % sort the axis
        k.ss.el.Y_s=k.ss.el.Y(k.ss.el.dist_idx);

        ss_id = bsxfun(@ge,k.ss.el.X_s,abs(k.qs2(:,1))') & bsxfun(@ge,k.ss.el.Y_s,abs(k.qs2(:,2))');
        k.ss.el.X_f=zeros(sum(ss_id(:,1)),4); 
        k.ss.el.Y_f=k.ss.el.X_f; 
        k.ss.el.dist_f=k.ss.el.X_f; 
        
        for q=1:4
            k.ss.el.X_f(:,q) = k.qs(q,1) * k.ss.el.X_s(ss_id(:,q)); 
            k.ss.el.Y_f(:,q) = k.qs(q,2) * k.ss.el.Y_s(ss_id(:,q));
            k.ss.el.dist_f(:,q) = k.ss.el.dist_s(ss_id(:,q));
        end
    end
    
    %% Generate the order for visting cells
    % Randomly permute the cell not known (to be visited). And generate the
    % Y.x and Y.y coordinate in a random order.
    
    if scale_i<parm.cstk_s
        parm.cstk=0;
    else
        parm.cstk=1;
    end
    
    Y{scale_i}.sim.xy=grid{scale_i}.xy(isnan(Y{scale_i}.m{1}));
    Y{scale_i}.sim.n=numel(Y{scale_i}.sim.xy);
    if parm.cstk
        Y{scale_i}.sim.xy_r=Y{scale_i}.sim.xy(randperm(Y{scale_i}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
        [Y{scale_i}.sim.x_r,Y{scale_i}.sim.y_r] = ind2sub([grid{scale_i}.ny, grid{scale_i}.nx],Y{scale_i}.sim.xy_r); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
    else
        for i_realisation=1:parm.n_realisation
            Y{scale_i}.sim.xy_r{i_realisation}=Y{scale_i}.sim.xy(randperm(Y{scale_i}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
            [Y{scale_i}.sim.x_r{i_realisation}, Y{scale_i}.sim.y_r{i_realisation}] = ind2sub([grid{scale_i}.ny, grid{scale_i}.nx],Y{scale_i}.sim.xy_r{i_realisation}); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
        end
    end

    %%
    % * *RANDOM NORMAL FIELD* 
    % This create the random Normal distribution used for sampling the posteriori distribution at each point.
    U = normcdf(randn(Y{scale_i}.sim.n,parm.n_realisation)); % random field



    %% Point simulation

    i_plot=0; % used for computing the number of point simulated. Used for ploting
    for i_pt=1:Y{scale_i}.sim.n; % loop over each point

        i_plot=i_plot+1;

        if parm.cstk % if option constant weight is activate.
            Y{scale_i}.pt.y = Y{scale_i}.sim.x_r(i_pt); % Find the current point position on the grid{scale_i}. It will be used on each realisation.
            Y{scale_i}.pt.x = Y{scale_i}.sim.y_r(i_pt);
            k0 = kringing_coef(Y{scale_i},X,k,parm,1); % Kriging.
            % variance of the kriging
            Y{scale_i}.pt.s = k0.s;
            
            % * *LIKELIHOOD*
            % We first find the secondary value (and error) (Z.d, Z.std) to
            % create a pdf. This pdf is then multiply inside the kernel to get
            % the density. The likelihood is only...
            Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x))  , Z.std(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x)));
            Y{scale_i}.pt.dens = bsxfun(@times, kernel.dens, Z.pt.dist./sum(Z.pt.dist));
            Y{scale_i}.pt.likelihood = sum(Y{scale_i}.pt.dens, 2)/sum(Y{scale_i}.pt.dens(:));
        end
        
       

        for i_realisation=1:parm.n_realisation

            %%
            % * *KRIGING*        
            if ~parm.cstk
                Y{scale_i}.pt.y = Y{scale_i}.sim.x_r{i_realisation}(i_pt); % Find the current point position on the grid{scale_i}. changing for each realisation
                Y{scale_i}.pt.x = Y{scale_i}.sim.y_r{i_realisation}(i_pt);
                k0 = kringing_coef(Y{scale_i},X,k,parm,i_realisation);
                % variance of the kriging
                Y{scale_i}.pt.s = k0.s;
                
                % * *LIKELIHOOD*
                % We first find the secondary value (and error) (Z.d, Z.std) to
                % create a pdf. This pdf is then multiply inside the kernel to get
                % the density. The likelihood is only...
                Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x))  , Z.std(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x)));
                Y{scale_i}.pt.dens = bsxfun(@times, kernel.dens, Z.pt.dist./sum(Z.pt.dist));
                Y{scale_i}.pt.likelihood = sum(Y{scale_i}.pt.dens, 2)/sum(Y{scale_i}.pt.dens(:));
            end

            if parm.neigh % if option smart neihbour is selected
                % the weight(lambda) were computed earlier in kriging_coef.m,
                % then the point are coming from the hard data (X.d_ns) and the
                % previously simulated point Y{scale_i}.m_ns.
                Y{scale_i}.pt.m = k0.lambda'* [X.d_ns(k0.sb_mask) ; Y{scale_i}.m_ns{i_realisation}(k0.ss_mask)];
            else
                XY_ns = [X.d_ns; Y{scale_i}.m_ns{i_realisation}(~isnan(Y{scale_i}.m_ns{i_realisation}))];
                Y{scale_i}.pt.m = k0.lambda'* XY_ns(k0.mask);
            end
            
            assert(~isnan(Y{scale_i}.pt.m),'The result can''t be NaN')


            %%
            % * * PRIOR DISTRIBUTION*
            % Back transform the normal distribution (from krigeage) in original space in the kernel.y grid{scale_i}. this become the prior distribution
            Y{scale_i}.pt.prior = Nscore.dist(Y{scale_i}.pt.m, sqrt(Y{scale_i}.pt.s));
            Y{scale_i}.pt.prior = Y{scale_i}.pt.prior./sum(Y{scale_i}.pt.prior);
            
%             x=X.d;
%             Nscore = nscore_perso(x, 'linear', 'linear', kernel, parm.plot.ns);
%             
%             figure; hold on;
%             [f,x2]=hist(X.d,kernel.y); plot(x2,f/sum(f),'linewidth',2);
%             plot(kernel.y,Nscore.dist(0,1))
%             
%             figure; hold on;
%             ecdf(Nscore.forward(X.d)); plot(Nscore.forward(kernel.y),cumsum(Nscore.dist(0,1)))
            
            
            %%
            % * *POSTERIORI*
            % Multiply the (nomalized) prior and likelihood to get the (normalised) posteriori
            Y{scale_i}.pt.post_pdf = Y{scale_i}.pt.prior.^(1-parm.p_w(scale_i)).* Y{scale_i}.pt.likelihood.^(parm.p_w(scale_i));
            Y{scale_i}.pt.post_pdf = Y{scale_i}.pt.post_pdf./sum(Y{scale_i}.pt.post_pdf);

            %%
            % * *SAMPLING*: Sample a point in the posteri distribution. To
            % do so the CDF is created and interpolation tool is used to assure
            % randomized sampling. 
            Y{scale_i}.pt.post_cdf = cumsum(Y{scale_i}.pt.post_pdf);

            % Interpolation only works if it is mono.. increasing. But the cdf
            % can reach 1 so we just add eps (very small value).
            if ~all(diff(Y{scale_i}.pt.post_cdf)>0) 
                Y{scale_i}.pt.post_cdf = Y{scale_i}.pt.post_cdf +  linspace(0,numel(Y{scale_i}.pt.post_cdf)*eps*2,numel(Y{scale_i}.pt.post_cdf))';
                assert(all(diff(Y{scale_i}.pt.post_cdf)>0))
            end

            Y{scale_i}.pt.sampled = interp1(Y{scale_i}.pt.post_cdf, kernel.y, U(i_pt,i_realisation),'pchip');


            %%
            % * *PLOTIT*
            if parm.plot.krig && i_realisation==1 && (i_plot==1 || mod(i_plot+49,50)==0  || Nscore.forward(Y{scale_i}.pt.sampled)<-3  || Nscore.forward(Y{scale_i}.pt.sampled)>3 ) % 
                figure(1); clf

                subplot(3,2,[1 4]);hold on
                h1=imagesc(Y{scale_i}.x,Y{scale_i}.y,Y{scale_i}.m_ns{1},'AlphaData',~isnan(Y{scale_i}.m_ns{1}));

                sb_i = min([round((Y{scale_i}.y(Y{scale_i}.pt.y)-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny]);
                sb_j = min([round((Y{scale_i}.x(Y{scale_i}.pt.x) -k.sb.x(1))/k.sb.dx +1)'; k.sb.nx]);
                windows=nan(k.sb.ny,k.sb.nx);
                for u=1:length(k.el_X_s) %.. look at all point...
                    for q=1:4 %... and assign it to the corresponding quadrant
                        if sb_i+k.qs(q,1)*k.el_X_s(u)<=k.sb.nx && sb_j+k.qs(q,2)*k.el_Y_s(u)<=k.sb.ny && sb_i+k.qs(q,1)*k.el_X_s(u)>=1 && sb_j+k.qs(q,2)*k.el_Y_s(u)>=1% check to be inside the grid
                            windows(sb_i+k.qs(q,1)*k.el_X_s(u), sb_j+k.qs(q,2)*k.el_Y_s(u))=1;
                        end
                    end
                end
                h2=imagesc(k.sb.x,k.sb.y,windows,'AlphaData',windows*.5);
                h3=mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.ny+1, k.sb.nx+1),'EdgeColor','k','facecolor','none');

                tt=-pi:0.01:pi;
                x=Y{scale_i}.x(Y{scale_i}.pt.x)+k.range(1)*cos(tt);
                x2=Y{scale_i}.x(Y{scale_i}.pt.x)+k.wradius*k.range(1)*cos(tt);
                y=Y{scale_i}.y(Y{scale_i}.pt.y)+k.range(2)*sin(tt);
                y2=Y{scale_i}.y(Y{scale_i}.pt.y)+k.wradius*k.range(2)*sin(tt);
                h4=plot(x,y,'--r'); h5=plot(x2,y2,'-r');
                h6=plot([Y{scale_i}.x(Y{scale_i}.pt.x) Y{scale_i}.x(Y{scale_i}.pt.x)],[min(y2) max(y2)],'-r');
                h7=plot([min(x2) max(x2)], [Y{scale_i}.y(Y{scale_i}.pt.y) Y{scale_i}.y(Y{scale_i}.pt.y)],'-r');

                lambda_c= 36+60.*(abs(k0.lambda)-min(abs(k0.lambda)))./range(abs(k0.lambda));

                h8=scatter(X.x,X.y,[],X.d_ns,'s','filled');

                if parm.cstk
                    n_hd = numel(X.x(k0.sb_mask));
                    sel_g=[X.x(k0.sb_mask) X.y(k0.sb_mask); Y{scale_i}.X(k0.ss_mask) Y{scale_i}.Y(k0.ss_mask)];
                    XY_ns = [X.d_ns(k0.sb_mask) ; Y{scale_i}.m_ns{i_realisation}(k0.ss_mask)];
                    h9=scatter(sel_g(1:n_hd,1),sel_g(1:n_hd,2),lambda_c(1:n_hd),XY_ns(1:n_hd),'s','filled','MarkerEdgeColor','k');
                    h10=scatter(sel_g(n_hd+1:end,1),sel_g(n_hd+1:end,2),lambda_c(n_hd+1:end),XY_ns(n_hd+1:end),'o','filled','MarkerEdgeColor','k');
                    h11=scatter(Y{scale_i}.x(Y{scale_i}.pt.x),Y{scale_i}.y(Y{scale_i}.pt.y),100,k0.lambda'* XY_ns,'o','filled','MarkerEdgeColor','r','LineWidth',1.5);
                else
                    XY_ns = [X.d_ns; Y{scale_i}.m_ns{1}(~isnan(Y{scale_i}.m_ns{1}))];
                    sel_g_ini=[X.x X.y; Y{scale_i}.X(~isnan(Y{scale_i}.m{i_realisation})) Y{scale_i}.Y(~isnan(Y{scale_i}.m{i_realisation}))];
                    sel_g = sel_g_ini(k0.mask,:);
                    scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns(k0.mask),'o','filled','MarkerEdgeColor','k');
                    scatter(Y{scale_i}.x(Y{scale_i}.pt.x),Y{scale_i}.y(Y{scale_i}.pt.y),100,k0.lambda'* XY_ns(k0.mask),'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
                end
                xlabel('x[m]');ylabel('y[m]');colorbar; 
                xlim([Y{scale_i}.x(1) Y{scale_i}.x(end)])
                ylim([Y{scale_i}.y(1) Y{scale_i}.y(end)])
                %set(gca,'YDir','reverse');
                
                legend([h3 h6 h8 h9 h10 h11],'Super grid','Window search with quadrant','Hard data point','Selected hard data point','Selected previously simulated point','Simulated Point','Location','northoutside','Orientation','horizontal')

                
                subplot(3,2,5); hold on;
                plot( kernel.y,Y{scale_i}.pt.prior)
                plot( kernel.y,Y{scale_i}.pt.likelihood)
                plot( kernel.y,Y{scale_i}.pt.post_pdf)
                plot(  kernel.y, max(Y{scale_i}.pt.post_pdf)*Y{scale_i}.pt.post_cdf)
                plot([Y{scale_i}.pt.sampled Y{scale_i}.pt.sampled],[0 max(Y{scale_i}.pt.post_pdf)],'k')
                legend(['Prior, w=' num2str(1-parm.p_w(scale_i))], ['Likelihood, w=' num2str(parm.p_w(scale_i))],'Posteriori','Post. cdf','sampled location')


                subplot(3,2,6); hold on;
                [f,x]=hist(X_ini.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                [f,x]=hist(Z.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                [f,x]=hist(Y{scale_i}.m{i_realisation}(:),50); plot(x,f/trapz(x,f));
                legend('Well sampling', 'ERT','Simulation')

                drawnow
                keyboard
            end

            %% Back transform 
            % * *
            Y{scale_i}.m{i_realisation}(Y{scale_i}.pt.y,Y{scale_i}.pt.x) = Y{scale_i}.pt.sampled;
            assert(~isnan(Nscore.forward(Y{scale_i}.pt.sampled)))
            Y{scale_i}.m_ns{i_realisation}(Y{scale_i}.pt.y,Y{scale_i}.pt.x) = Nscore.forward(Y{scale_i}.pt.sampled); % add only the last simulated point to normal score

        end
    end

   
    % display info
    disp(['Simulation ' num2str(scale_i) '/' num2str(parm.n_scale) ' finished in : ' num2str(toc(t.tic.scale))])
    t.scale{scale_i} = toc(t.tic.scale);
end

% save intial value
X = X_ini;
if parm.neigh; k.sb.mask = k.sb.mask_ini; end
clear X_ini k.sb.mask_ini

t.global = toc(t.tic.global);
clear t.tic


%% * 5. *SAVE IT*
if parm.saveit
    mkdir(['result/', parm.familyname])
    filename=['result/', parm.familyname, 'SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    save(filename, 'parm', 'Y', 'grid', 't', 'X', 'Z', 'X_true', 'k', 'kernel', 'Nscore')
end


end