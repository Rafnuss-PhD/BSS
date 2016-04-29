%% BSGS is an implementation of the Bayesian Sequential Gaussian Simulation
% BSGS is a stochatic sequential simulation of a primary variable (X) based
% on (1) some initialy-known scattered primary variable (hard data) and (2)
% a coarse secondary variable information (Z).
%
%
%
% INPUT:
%
% * Prim.d      : Primary variable (1D-X.n elts)
% * Prim.x,y    : x,y-coordinate of the primary variable (1D-X.n elts)
% * Sec.d       : Secondary variable (2D- grid.nx x grid.n elts)
% * Sec.x,y     : x,y-coordinate of the secondary variable (1D elts)
% * Sec.std     : Secondary variable std error (2D-grid.nx x grid.n elts)
% * grid_gen    : grid informations
% * parm        : Parameters of the simulation.
%   * gen       : generation parmater structure (ouput of data_generation.m)
%   * likelihood: bolean. with or without using likelihood.
%   * scale     : array of the scales (i, where nx=2^grid{i}+1) to simulate
%   * name      : name of the simulation
%   * seed      : random number. use to reproduce exactly the same simulation.
%   * saveit 	: save or not in a file
%   * unit 	: unit of the primary variable, used for plot
%   * n_realisation : number of realisations.
%   * neigh 	: bolean. with or without "smart neighbour"
%   * nscore	: bolean. with or without  normal score transform
%   * cstk      : bolean. with or without constant path
%   * cstk_s 	: scale at which cst path is switch on
%   * plot:
%       * bsgs	: plot. ?
%       * ns 	: plot. ?
%       * sb 	: plot. ?
%       * kernel : plot. ?
%       * fitvar : plot. ?
%       * krig  : plot. ?
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
% * R       : Resultant Primary variable
% * t       : Time of simulations
% * kernel  : kernel information
% * k       : kriging information
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% Referances:
%

function [Res, t, k, parm, filename] = SGSIM(Prim,Prim_true,grid_gen,parm)
t.global = tic;
addpath(genpath('./.'))

%% * *INPUT CEHCKING*
% This section of the code generates a valid parm structure based on the
% inputted parm. If some parm value haven't been input, this section will
% fill automatically with defautl value. This may not allowed be the best.

% Force input in column
Prim.x=Prim.x(:);Prim.y=Prim.y(:);Prim.d=Prim.d(:);

% Parameter settings
if ~isfield(parm, 'seed'),          parm.seed           = rand(); end
rng(parm.seed);
if ~isfield(parm, 'saveit'),        parm.saveit         = 1; end % bolean, save or not the result of simulation
if ~isfield(parm, 'name'),          parm.name           = ''; end % name use for saving file
if ~isfield(parm, 'familyname'),    parm.familyname     = ''; end
if ~isfield(parm, 'unit'),          parm.unit           = ''; end % unit of the primary variable, used in plot
if ~isfield(parm, 'n_realisation'), parm.n_realisation  = 1; end
if ~isfield(parm, 'par'),           parm.par            = 1; end
if ~isfield(parm, 'par_n'),         parm.par_n          = feature('numcores'); end
if ~isfield(parm, 'notify'),
    parm.notify          = 0;
else
    if ~isfield(parm, 'notify_email'), parm.notify_email  = 'rafnuss@gmail.com'; end
end
if ~isfield(parm, 'path'),         parm.path            = 'random'; end
if ~isfield(parm, 'varcovar'),         parm.path            = 0; end
% Scale and weight parameters
if ~isfield(parm, 'scale')
    parm.scale = repmat(1:max([grid_gen.sx,grid_gen.sy]),2,1);
    parm.scale(1,parm.scale(1,:)>grid_gen.sx) = grid_gen.sx;
    parm.scale(2,parm.scale(2,:)>grid_gen.sy) = grid_gen.sy;
end
if ~isfield(parm, 'cstk_s') % cstk_s is the scale at which cst is switch on
    if ~isfield(parm, 'cstk'),      parm.cstk           = 1; end % constant path and kriging weight activated or not
    if parm.cstk
        parm.cstk_s = 0; % will never use cstk
    else
        parm.cstk_s = Inf; % will always use cstk
    end
end
if ~isfield(parm, 'fitvar'), parm.fitvar     =0; end % fit the variogram to the data or used the given one in parm.covar

% Nscore
parm.support_dist = linspace(-5,5,500)';

% Kriging parameter
if ~isfield(parm, 'neigh'),         parm.neigh          =1; end % smart-neighbouring activated or not
if ~isfield(parm, 'nscore'),        parm.nscore         =1; end % use normal score (strongly advice to use it.)
if ~isfield(parm, 'k') || ~isfield(parm.k, 'nb_neigh')
    parm.k.nb_neigh = [0 0 0 0 0; 5 5 5 5 5];
end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'model') || ~isfield(parm.k, 'c')
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
k = parm.k;


% Plot
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'bsgs'),  parm.plot.bsgs   =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'ns'),    parm.plot.ns     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'sb'),    parm.plot.sb     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'kernel'),parm.plot.kernel =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'fitvar'),parm.plot.fitvar =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'krig'),  parm.plot.krig   =0; end


% Compute the rot matrix
for i=1:size(parm.k.model,1)
    ang=parm.k.model(i,4); cang=cos(ang/180*pi); sang=sin(ang/180*pi);
    rot = [cang,-sang;sang,cang];
    parm.k.cx{i} = rot/diag(parm.k.model(i,2:3));
end

% Check the input for correct size dans dimension
assert(size(Prim.d,2)<=1,'X.d is not a vertical vector 1D');
assert(size(Prim.x,2)<=1,'X.dx is not a vertical vector 1D');
assert(size(Prim.y,2)<=1,'X.y is not a vertical vector 1D');
assert(all(size(Prim.y)==size(Prim.x)),'X.x and X.y don''t have the same dimension');
assert(all(size(Prim.d)==size(Prim.x)),'X.d and X.x (or X.y) don''t have the same dimension');
assert(max(parm.scale(:))<=max([grid_gen.sx,grid_gen.sy]),'This scale of simulation does not exist')
assert(all(parm.scale(1,:)<=grid_gen.sx),'monotonicly increasing scale')
assert(all(parm.scale(2,:)<=grid_gen.sy),'monotonicly increasing scale')

% Creation of the grid
parm.n_scale=size(parm.scale,2);
grid=cell(parm.n_scale,1);
for i_scale = 1:parm.n_scale
    grid{i_scale}.sx=parm.scale(1,i_scale);
    grid{i_scale}.sy=parm.scale(2,i_scale);
    grid{i_scale}.nx=2^grid{i_scale}.sx+1;
    grid{i_scale}.ny=2^grid{i_scale}.sy+1;
    grid{i_scale}.nxy=grid{i_scale}.nx*grid{i_scale}.ny; % total number of cells
    
    grid{i_scale}.dx=grid_gen.x(end)/(grid{i_scale}.nx-1);
    grid{i_scale}.dy=grid_gen.y(end)/(grid{i_scale}.ny-1);
    
    grid{i_scale}.x=linspace(0, grid_gen.x(end), grid{i_scale}.nx); % coordinate of cells center
    grid{i_scale}.y=linspace(0, grid_gen.y(end), grid{i_scale}.ny);
    grid{i_scale}.xy=1:grid{i_scale}.nxy;
    
    [grid{i_scale}.X, grid{i_scale}.Y] = meshgrid(grid{i_scale}.x,grid{i_scale}.y); % matrix coordinate
end


%% * 1. *SUPERBLOCK GRID CREATION*
% A mask (Boolean value) of the hard data is assigned to each superblock
% as follow: Only the n-closest (normalized by the covariance range) points
% (inside the ellipse/windows) to the centre of the superblock will be
% true. During the kriging, the mask of the superblock of the estimated
% point will be used to select the hard to add to the kriging system
if parm.neigh
    k.sb.nx = parm.k.sb.nx; % number of superblock grid
    k.sb.ny = parm.k.sb.ny;
    [k, Prim] = SuperBlockGridCreation(k, grid_gen.x(end), grid_gen.y(end), Prim, parm.plot.sb, parm.k.nb_neigh(2,5));
end


% %% * 2. *NON-PARAMETRIC RELATIONSHIP*
% % The joint pdf of the primary and secondary is build using a bivariate
% % kernel density estimator (Botev, Grotowski, & Kroese, 2010).
% if Prim.n~=0
%     kern = kernel(Prim, Sec, parm.kernel_range, parm.plot.kernel);
% else
%     warning('No hard data !')
%     assert(all(parm.p_w(1,:)==0),'The Secondary cannot be used !')
%     kern.axis_prim=linspace(-5,5,100)';
%     kern.axis_sec=kern.axis_prim;
%     kern.dens=ones(100,100);
% end


%% * 3. *NORMAL SCORE TRANSFORM*
% Based on the hard data (well samples), a normal score transform is
% created using interpolation with power or exponential tail extrapolation.
% Matlab symbolique function are used for efficient coding. The back
% transform of the prior normal distribution function is also created
% (return the pdf in the initial space from the mean and variance in the
% normal space)
if parm.nscore && Prim.n~=0
    Nscore = nscore(Prim.d, parm.support_dist, 'linear', 'linear', parm.plot.ns);
else
    Nscore.forward = @(x) x;
    Nscore.inverse = @(x) x;
    Nscore.dist    = @(mu,sigma) normpdf(parm.support_dist,mu,sigma)/sum(normpdf(parm.support_dist,mu,sigma));
end

% Create the normal space primary variable of known data
Prim.d_ns = Nscore.forward(Prim.d);



%% * 4. *RUN SIMULATION*

if parm.neigh; k.sb.mask_ini = k.sb.mask; end % mask will be change after, we preserve its structure

if parm.par && parm.n_realisation~=1 % if parralelelisation is selected
    delete(gcp('nocreate'));
    poolobj=parpool(parm.par_n); % find the number of core available
    par_n_realisation = ceil(parm.n_realisation/poolobj.NumWorkers);
    
    RR=cell(poolobj.NumWorkers,1);
    tt=cell(poolobj.NumWorkers,1);
    
    parm_pool=parm;
    parm_pool.n_realisation = par_n_realisation;
    
    parfor pool_i=1:poolobj.NumWorkers
        [RR{pool_i}, tt{pool_i}]=SGSIM_in(Prim, k, Nscore, grid, parm_pool);
    end
    delete(poolobj)
    
    Res=RR{1};
    for pool_i=2:numel(RR)
        for i_scale=1:parm.n_scale
            Res{i_scale}.m = [Res{i_scale}.m; RR{pool_i}{i_scale}.m];
            Res{i_scale}.m_ns = [Res{i_scale}.m_ns; RR{pool_i}{i_scale}.m_ns];
        end
    end
    
    
else
    [Res,t.scale,t.cstk] = SGSIM_in(Prim, k, Nscore, grid, parm);
end


% save intial value
if parm.neigh; k.sb.mask = k.sb.mask_ini; end
clear X_ini k.sb.mask_ini


t.global = toc(t.global );


%% * 5. *SAVE IT*
filename=['result-SGSIM/', parm.familyname, 'SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
if parm.saveit
    mkdir(['result-SGSIM/', parm.familyname])
    save(filename, 'parm', 'Res', 'grid', 't', 'Prim', 'Prim_true', 'k',  'Nscore')
end

if parm.notify
    unix(['echo "Simulation ' parm.familyname ' has finish now (' datestr(now,'yyyy-mm-dd_HH-MM-SS') ') in ' num2str(t.global/60) 'min" | mail -s "Simulation Finish" ' parm.notify_email]);
    load handel; sound(y,Fs)
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Res,t_scale,t_cstk]=SGSIM_in(Prim, krig, Nscore, grid, parm)

Res = cell(parm.n_scale,1);
t_scale = zeros(parm.n_scale,1);
t_cstk = zeros(parm.n_scale,1);

if parm.varcovar % compute the empirical var-covariance matrix (Emery & Pelaez, 2011)
    Prim.varcovar=zeros(Prim.n,1);
    
    % Assign to each hard data the index of the final grid
    for i_hard_data=1:Prim.n
        if ismember(Prim.y(i_hard_data),grid{parm.n_scale}.y) && ismember(Prim.x(i_hard_data),grid{parm.n_scale}.x) % if belong to the grid
            Prim.varcovar(i_hard_data) = find(grid{parm.n_scale}.Y==Prim.y(i_hard_data) & grid{parm.n_scale}.X==Prim.x(i_hard_data));
        end
    end
    
    Prim.varcovar_safe= Prim.varcovar;
    Res{1}.lambda0 = sparse(grid{parm.n_scale}.nxy, sum(Prim.varcovar==0)); % for point not on the grid
    Res{1}.lambda = sparse(grid{parm.n_scale}.nxy,grid{parm.n_scale}.nxy);
    
end

%% * 1. *Simulation of Scale*
for i_scale=1:parm.n_scale % for each scale
    t.tic.scale = tic;
    
    %% * 1.1. *INITIATE SCALE SIMULATION*
    % Allocating space for resulting field.
    Res{i_scale}.x=grid{i_scale}.x;
    Res{i_scale}.y=grid{i_scale}.y;
    Res{i_scale}.X=grid{i_scale}.X;
    Res{i_scale}.Y=grid{i_scale}.Y;
    Res{i_scale}.nx=grid{i_scale}.nx;
    Res{i_scale}.ny=grid{i_scale}.ny;
    Res{i_scale}.nxy=grid{i_scale}.ny*Res{i_scale}.nx;
    Res{i_scale}.m=repmat({nan(grid{i_scale}.ny,grid{i_scale}.nx)},parm.n_realisation,1); % matrix des resutlats
    
    % Populate the grid from previous scale.
    if i_scale~=1 % not at first scale
        for i_sim=1:parm.n_realisation
            Res{i_scale}.m{i_sim}( 1:(grid{i_scale-1}.dy/grid{i_scale}.dy):end, 1:(grid{i_scale-1}.dx/grid{i_scale}.dx):end) = Res{i_scale-1}.m{i_sim};
        end
    end
    
    % Assimilate the hard data (Prim) into the grid
    hard_data_idx=find(ismember(Prim.y,grid{i_scale}.Y)&ismember(Prim.x,grid{i_scale}.X));
    for i_sim=1:parm.n_realisation
        for hd=1:numel(hard_data_idx)
            Res{i_scale}.m{i_sim}(Prim.x(hard_data_idx(hd))==grid{i_scale}.X & Prim.y(hard_data_idx(hd))==grid{i_scale}.Y) = Prim.d(hard_data_idx(hd));
        end
    end
    
    % Remove the assimilated data from X.
    Prim.d(hard_data_idx)=[];
    Prim.x(hard_data_idx)=[];
    Prim.y(hard_data_idx)=[];
    Prim.d_ns(hard_data_idx)=[];
    if parm.varcovar; Prim.varcovar(hard_data_idx)=[]; end
    if parm.neigh; krig.sb.mask(:,:,hard_data_idx)=[]; end
    Prim.n=numel(Prim.d);
    
    % Create the normal space result matrix and populate with known value
    Res{i_scale}.m_ns = repmat({nan(size(Res{i_scale}.m{1}))},parm.n_realisation,1);
    I=~isnan(Res{i_scale}.m{1});
    for i_realisation=1:parm.n_realisation
        Res{i_scale}.m_ns{i_realisation}(I) = Nscore.forward(Res{i_scale}.m{i_realisation}(I));
    end
    

    if parm.varcovar
        Res{i_scale}.varcovar_id=reshape( find(ismember(grid{parm.n_scale}.X, Res{i_scale}.x) &  ismember(grid{parm.n_scale}.Y, Res{i_scale}.y)), Res{i_scale}.ny, Res{i_scale}.nx);
    end
    
    
    
    %% * 1.2. *SPIRAL SEARCH*
    % Create the windows for kringing with the function to compute the
    % normalized distence and the order of visit of the cells. Spiral
    % Search setting: previously data (on grid{i_scale} location)
    if parm.neigh
        [krig.ss.el.X, krig.ss.el.Y] = meshgrid(0:max(ceil(krig.range(1)*krig.wradius/grid{i_scale}.dx),ceil(krig.range(2)*krig.wradius/grid{i_scale}.dy)));% grid{i_scale} of searching windows
        [krig.ss.el.X_T, krig.ss.el.Y_T]=rotredtrans(krig.ss.el.X*grid{i_scale}.dx, krig.ss.el.Y*grid{i_scale}.dy, krig.rotation, krig.range); % transforms the grid{i_scale}
        krig.ss.el.dist = sqrt(krig.ss.el.X_T.^2 + krig.ss.el.Y_T.^2); % find distence
        [krig.ss.el.dist_s, krig.ss.el.dist_idx] = sort(krig.ss.el.dist(:)); % sort according distence.
        krig.ss.el.X_s=krig.ss.el.X(krig.ss.el.dist_idx); % sort the axis
        krig.ss.el.Y_s=krig.ss.el.Y(krig.ss.el.dist_idx);
        
        ss_id = bsxfun(@ge,krig.ss.el.X_s,abs(krig.qs2(:,1))') & bsxfun(@ge,krig.ss.el.Y_s,abs(krig.qs2(:,2))');
        krig.ss.el.X_f=zeros(sum(ss_id(:,1)),4);
        krig.ss.el.Y_f=krig.ss.el.X_f;
        krig.ss.el.dist_f=krig.ss.el.X_f;
        
        for q=1:4
            krig.ss.el.X_f(:,q) = krig.qs(q,1) * krig.ss.el.X_s(ss_id(:,q));
            krig.ss.el.Y_f(:,q) = krig.qs(q,2) * krig.ss.el.Y_s(ss_id(:,q));
            krig.ss.el.dist_f(:,q) = krig.ss.el.dist_s(ss_id(:,q));
        end
    end
    
    
    %% * 1.3. *Generate the order for visting cells*
    % Randomly permute the cell not known (to be visited). And generate the
    % Y.x and Y.y coordinate in a random order.
    
    if i_scale<parm.cstk_s
        parm.cstk=0;
    else
        parm.cstk=1;
    end
    
    Res{i_scale}.sim.xy=grid{i_scale}.xy(isnan(Res{i_scale}.m{1}));
    Res{i_scale}.sim.n=numel(Res{i_scale}.sim.xy);
    if parm.cstk
        switch parm.path
            case 'random'
                Res{i_scale}.path = randperm(Res{i_scale}.sim.n);
            case 'row-by-row'
                Res{i_scale}.path = 1:Res{i_scale}.sim.n;
            case 'spiralin'
                dist = Inf*ones(Res{i_scale}.sim.n,1);
                X = Res{i_scale}.X(~isnan(Res{i_scale}.m{1}));
                Y = Res{i_scale}.Y(~isnan(Res{i_scale}.m{1}));
                for i_hard_data=1:numel(X)
                    dist_chall = sqrt( (Res{i_scale}.X(Res{i_scale}.sim.xy)-X(i_hard_data)).^2 + (Res{i_scale}.Y(Res{i_scale}.sim.xy)-Y(i_hard_data)).^2 );
                    dist = min(dist,dist_chall');
                end
                [~,Res{i_scale}.path] = sort(dist);
                clear dist dist_chall X Y
            case 'spiralout'
                dist = Inf*ones(Res{i_scale}.sim.n,1);
                X = Res{i_scale}.X(~isnan(Res{i_scale}.m{1}));
                Y = Res{i_scale}.Y(~isnan(Res{i_scale}.m{1}));
                for i_hard_data=1:numel(X)
                    dist_chall = sqrt( (Res{i_scale}.X(Res{i_scale}.sim.xy)-X(i_hard_data)).^2 + (Res{i_scale}.Y(Res{i_scale}.sim.xy)-Y(i_hard_data)).^2 );
                    dist = min(dist,dist_chall');
                end
                [~,Res{i_scale}.path] = sort(dist,'descend');
                clear dist dist_chall X Y
            case 'maximize'
                dist = Inf*ones(Res{i_scale}.sim.n,1);
                X = Res{i_scale}.X(~isnan(Res{i_scale}.m{1}));
                Y = Res{i_scale}.Y(~isnan(Res{i_scale}.m{1}));
                for i_hard_data=1:numel(X)
                    dist_chall = sqrt( (Res{i_scale}.X(Res{i_scale}.sim.xy)-X(i_hard_data)).^2 + (Res{i_scale}.Y(Res{i_scale}.sim.xy)-Y(i_hard_data)).^2 );
                    dist = min(dist,dist_chall');
                end
                Res{i_scale}.path = nan(Res{i_scale}.sim.n,1);
                for i_pt = 1:Res{i_scale}.sim.n
                    [~,Res{i_scale}.path(i_pt)] = max(dist);
                    dist_chall = sqrt( (Res{i_scale}.X(Res{i_scale}.sim.xy)-Res{i_scale}.X(Res{i_scale}.sim.xy(Res{i_scale}.path(i_pt)))).^2 + (Res{i_scale}.Y(Res{i_scale}.sim.xy)-Res{i_scale}.Y(Res{i_scale}.sim.xy(Res{i_scale}.path(i_pt)))).^2 );
                    dist = min(dist,dist_chall');
                end
                clear dist dist_chall X Y
            otherwise
                error('path method not define')
              
        end
        Res{i_scale}.sim.xy_r=Res{i_scale}.sim.xy(Res{i_scale}.path); % randomly permutate the ordered vector of index of Y.xy
        [Res{i_scale}.sim.y_r,Res{i_scale}.sim.x_r] = ind2sub([grid{i_scale}.ny, grid{i_scale}.nx],Res{i_scale}.sim.xy_r); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
        
    else
        for i_realisation=1:parm.n_realisation
            Res{i_scale}.sim.xy_r{i_realisation}=Res{i_scale}.sim.xy(randperm(Res{i_scale}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
            [Res{i_scale}.sim.y_r{i_realisation}, Res{i_scale}.sim.x_r{i_realisation}] = ind2sub([grid{i_scale}.ny, grid{i_scale}.nx],Res{i_scale}.sim.xy_r{i_realisation}); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
        end
    end
    
    
    %% * 1.4. *RANDOM NORMAL FIELD*
    % This create the random Normal distribution used for sampling the
    % posteriori distribution at each point.
    U = normcdf(randn(Res{i_scale}.sim.n,parm.n_realisation)); % random field
    
    
    %% * 2. *Simulation of Point*
    i_plot=0; % used for computing the number of point simulated. Used for ploting
    for i_pt=1:Res{i_scale}.sim.n; % loop over each point
        i_plot=i_plot+1;
        
%         if i_pt==492
%             deo
%         end
            
        if parm.cstk % if option constant weight is activate.
            % Find the current point position on the grid{i_scale}. It
            % will be used on each realisation.
            t.tic.cstk = tic;
            pt.x = Res{i_scale}.sim.x_r(i_pt);
            pt.y = Res{i_scale}.sim.y_r(i_pt);
            
            % Kriging system
            pt.krig = kriging(pt,Res{i_scale},Prim,krig,parm,1);
            t_cstk(i_scale)=t_cstk(i_scale)+toc(t.tic.cstk);
        end
        
        
        %% * 3. *Simulation of Realisation*
        for i_realisation=1:parm.n_realisation
            if ~parm.cstk
                % Find the current point position on the grid{i_scale}.
                % This is changing for each realisation
                pt.x = Res{i_scale}.sim.x_r{i_realisation}(i_pt);
                pt.y = Res{i_scale}.sim.y_r{i_realisation}(i_pt);
                
                % Kriging
                pt.krig = kriging(pt,Res{i_scale},Prim,krig,parm,i_realisation);
            end
            
            if parm.neigh % if option smart neihbour is selected
                % the weight(lambda) were computed earlier in kriging_coef.m,
                % then the point are coming from the hard data (X.d_ns) and the
                % previously simulated point R{i_scale}.m_ns.
                pt.krig.m = pt.krig.lambda'* [Prim.d_ns(pt.krig.sb_mask) ; Res{i_scale}.m_ns{i_realisation}(pt.krig.ss_mask)];
            else
                XY_ns = [Prim.d_ns; Res{i_scale}.m_ns{i_realisation}(~isnan(Res{i_scale}.m_ns{i_realisation}))];
                pt.krig.m = pt.krig.lambda'* XY_ns(pt.krig.mask);
            end
            
            if isempty(pt.krig.m)
                warning('Empty estimate. replace by 0')
                pt.krig.m=0;
            end
            
            %% * 3.1 *KRIGING DISTRIBUTION*
            % Back transform the normal distribution (from krigeage) in
            % original space in the  grid{i_scale}. this become
            % the prior distribution
            pt.krig.pdf = Nscore.dist(pt.krig.m, sqrt(pt.krig.s));
            pt.krig.pdf = pt.krig.pdf./sum(pt.krig.pdf);
            
            
            
            %% * 3.3  *SAMPLING*:
            % Sample a point in the posteri distribution. To do so the CDF
            % is created and interpolation tool is used to assure
            % randomized sampling.
            pt.krig.cdf = cumsum(pt.krig.pdf);
            
            % Interpolation only works if it is mono.. increasing. But the
            % cdf can reach 1 so we just add eps (very small value).
            if ~all(diff(pt.krig.cdf)>0)
                pt.krig.cdf = pt.krig.cdf +  linspace(0,numel(pt.krig.cdf)*eps*2,numel(pt.krig.cdf))';
                assert(all(diff(pt.krig.cdf)>0))
            end
            
            pt.sampled = interp1(pt.krig.cdf, parm.support_dist, U(i_pt,i_realisation),'pchip');
            
            
            %% * 3.4 *PLOTIT*
            if parm.plot.krig && i_realisation==1 && (i_plot==1261   || Nscore.forward(pt.sampled)<-3  || Nscore.forward(pt.sampled)>3 ) %|| mod(i_plot,50)==0
                figure(1); clf
                
                subplot(3,2,[1 4]);
                hold on
                h1=imagesc(Res{i_scale}.x,Res{i_scale}.y,Res{i_scale}.m_ns{1},'AlphaData',~isnan(Res{i_scale}.m_ns{1}));
                
                if parm.neigh
                    sb_i = min([round((Res{i_scale}.y(pt.y)-krig.sb.y(1))/krig.sb.dy +1)'; krig.sb.ny]);
                    sb_j = min([round((Res{i_scale}.x(pt.x) -krig.sb.x(1))/krig.sb.dx +1)'; krig.sb.nx]);
                    windows=nan(krig.sb.ny,krig.sb.nx);
                    for u=1:length(krig.el_X_s) %.. look at all point...
                        for q=1:4 %... and assign it to the corresponding quadrant
                            if sb_i+krig.qs(q,1)*krig.el_X_s(u)<=krig.sb.nx && sb_j+krig.qs(q,2)*krig.el_Y_s(u)<=krig.sb.ny && sb_i+krig.qs(q,1)*krig.el_X_s(u)>=1 && sb_j+krig.qs(q,2)*krig.el_Y_s(u)>=1% check to be inside the grid
                                windows(sb_i+krig.qs(q,1)*krig.el_X_s(u), sb_j+krig.qs(q,2)*krig.el_Y_s(u))=1;
                            end
                        end
                    end
                    h2=imagesc(krig.sb.x,krig.sb.y,windows,'AlphaData',windows*.5);
                    h3=mesh([0 krig.sb.x+krig.sb.dx/2],[0 krig.sb.y+krig.sb.dy/2],zeros(krig.sb.ny+1, krig.sb.nx+1),'EdgeColor','k','facecolor','none');
                end
                
                tt=-pi:0.01:pi;
                x=Res{i_scale}.x(pt.x)+krig.range(1)*cos(tt);
                x2=Res{i_scale}.x(pt.x)+krig.wradius*krig.range(1)*cos(tt);
                y=Res{i_scale}.y(pt.y)+krig.range(2)*sin(tt);
                y2=Res{i_scale}.y(pt.y)+krig.wradius*krig.range(2)*sin(tt);
                h4=plot(x,y,'--r'); h5=plot(x2,y2,'-r');
                h6=plot([Res{i_scale}.x(pt.x) Res{i_scale}.x(pt.x)],[min(y2) max(y2)],'-r');
                h7=plot([min(x2) max(x2)], [Res{i_scale}.y(pt.y) Res{i_scale}.y(pt.y)],'-r');
                
                lambda_c= 36+60.*(abs(pt.krig.lambda)-min(abs(pt.krig.lambda)))./range(abs(pt.krig.lambda));
                
                h8=scatter(Prim.x,Prim.y,[],Prim.d_ns,'s','filled');
                
                if parm.neigh
                    n_hd = numel(Prim.x(pt.krig.sb_mask));
                    sel_g=[Prim.x(pt.krig.sb_mask) Prim.y(pt.krig.sb_mask); Res{i_scale}.X(pt.krig.ss_mask) Res{i_scale}.Y(pt.krig.ss_mask)];
                    XY_ns = [Prim.d_ns(pt.krig.sb_mask) ; Res{i_scale}.m_ns{i_realisation}(pt.krig.ss_mask)];
                    h9=scatter(sel_g(1:n_hd,1),sel_g(1:n_hd,2),lambda_c(1:n_hd),XY_ns(1:n_hd),'s','filled','MarkerEdgeColor','k');
                    h10=scatter(sel_g(n_hd+1:end,1),sel_g(n_hd+1:end,2),lambda_c(n_hd+1:end),XY_ns(n_hd+1:end),'o','filled','MarkerEdgeColor','k');
                    h11=scatter(Res{i_scale}.x(pt.x),Res{i_scale}.y(pt.y),100,pt.krig.lambda'* XY_ns,'o','filled','MarkerEdgeColor','r','LineWidth',1.5);
                else
                    XY_ns = [Prim.d_ns; Res{i_scale}.m_ns{1}(~isnan(Res{i_scale}.m_ns{1}))];
                    sel_g_ini=[Prim.x Prim.y; Res{i_scale}.X(~isnan(Res{i_scale}.m{i_realisation})) Res{i_scale}.Y(~isnan(Res{i_scale}.m{i_realisation}))];
                    sel_g = sel_g_ini(pt.krig.mask,:);
                    scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns(pt.krig.mask),'o','filled','MarkerEdgeColor','k');
                    scatter(Res{i_scale}.x(pt.x),Res{i_scale}.y(pt.y),100,pt.krig.lambda'* XY_ns(pt.krig.mask),'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
                end
                xlabel('x[m]');ylabel('y[m]');
                %colorbar;
                xlim([Res{i_scale}.x(1) Res{i_scale}.x(end)])
                ylim([Res{i_scale}.y(1) Res{i_scale}.y(end)])
                %set(gca,'YDir','reverse');
                
                % legend([h3 h6 h8 h9 h10 h11],'Super grid','Window search with quadrant','Hard data point','Selected hard data point','Selected previously simulated point','Simulated Point','Location','northoutside','Orientation','horizontal')
                
                
                subplot(3,2,5); hold on;
                plot( parm.support_dist,pt.krig.pdf)
                plot(  parm.support_dist, max(pt.krig.pdf)*pt.krig.cdf)
                plot([pt.sampled pt.sampled],[0 max(pt.krig.pdf)],'k')
                legend('Kiriging estimate','CDF','sampled')
                
                
                subplot(3,2,6); hold on;
                %[f,x]=hist(Prim_ini.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                %[f,x]=hist(Z.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                [f,x]=hist(R{i_scale}.m{i_realisation}(:),50); plot(x,f/trapz(x,f));
                legend('Well sampling', 'ERT','Simulation')
                
                drawnow
                keyboard
            end
            
                    
            %% 3.5 Compute the empirical var-covariance matrix (Emery & Pelaez, 2011)
            if parm.varcovar && i_realisation==1
                pt_id = find(Res{i_scale}.y(pt.y)==grid{parm.n_scale}.Y&Res{i_scale}.x(pt.x)==grid{parm.n_scale}.X);
                
                if parm.neigh % if option smart neihbour is selected
                    Res{end}.lambda(pt_id, [Prim.varcovar(pt.krig.sb_mask); Res{i_scale}.varcovar_id(pt.krig.ss_mask)]) = -pt.krig.lambda./sqrt(pt.krig.s);
                    Res{end}.lambda(pt_id,pt_id) = 1/sqrt(pt.krig.s);
                else % compute only the last realization
                    XY_id = [Prim.varcovar; Res{i_scale}.varcovar_id(~isnan(Res{i_scale}.m{i_realisation}))];
                    Res{end}.lambda(pt_id,XY_id(pt.krig.mask)) = -pt.krig.lambda./sqrt(pt.krig.s);
                    Res{end}.lambda(pt_id,pt_id) = 1/sqrt(pt.krig.s);
                    clear XY_id
                end
                clear pt_id
            end
            
            %% 3.6 *Back transform*
            Res{i_scale}.m{i_realisation}(pt.y,pt.x) = pt.sampled;
            Res{i_scale}.m_ns{i_realisation}(pt.y,pt.x) = Nscore.forward(pt.sampled); % add only the last simulated point to normal score
        end
        % NOTHING SHOULD BE DONE AFTER THAT WHICH INVOLVED THE WEIGHT
        % BECAUSE WE CHANGED RES.M...
    end
    
    % display info
    disp(['Simulation ' num2str(i_scale) '/' num2str(parm.n_scale) ' finished in : ' num2str(toc(t.tic.scale))])
    t_scale(i_scale) = toc(t.tic.scale);
    

end

end
