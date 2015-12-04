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

function [Y, t, kernel, k] = BSGS_par(X,Z,X_true,grid_gen,parm)
t.global = tic;
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
if ~isfield(parm, 'par'),           parm.par            = 1; end
if ~isfield(parm, 'par_n'),         parm.par_n          = feature('numcores'); end
if ~isfield(parm, 'scale')
    parm.scale = repmat(1:max([grid_gen.sx,grid_gen.sy]),2,1); 
    parm.scale(1,parm.scale(1,:)>grid_gen.sx) = grid_gen.sx;
    parm.scale(2,parm.scale(2,:)>grid_gen.sy) = grid_gen.sy;
end
if ~isfield(parm, 'p_w')
    parm.p_w     = repmat(0.5, 2, size(parm.scale,2)); 
elseif size(parm.p_w,1)==1
    parm.p_w     = [parm.p_w ; 1-parm.p_w ];
end
if size(parm.p_w,2)==1
    parm.p_w     = repmat(parm.p_w, 1,size(parm.scale,2)); 
elseif size(parm.p_w,2)==2
    parm.p_w     = [linspace(parm.p_w(1,1), parm.p_w(1,2),size(parm.scale,2)); linspace(parm.p_w(2,1), parm.p_w(2,2),size(parm.scale,2))];
end
assert(size(parm.p_w,1)==2 & size(parm.p_w,2)==size(parm.scale,2))

% Run option
if ~isfield(parm, 'neigh'),         parm.neigh          =1; end % smart-neighbouring activated or not
if ~isfield(parm, 'nscore'),        parm.nscore         =1; end % use normal score (strongly advice to use it.)
if ~isfield(parm, 'cstk_s')% cstk_s is the scale at which cst is switch on
    if ~isfield(parm, 'cstk'),      parm.cstk           = 1; end % constant path and kriging weight activated or not
    if parm.cstk
        parm.cstk_s = 0; % will never use cstk
    else
        parm.cstk_s = Inf; % will always use cstk
    end
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

if parm.fitvar
    [k.range, k.var] = fit_variogramm(X,Z,parm.plot.fitvar, X_true);
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

% Create the normal space primary variable of known data
X.d_ns = Nscore.forward(X.d);



%% * 4. *SIMULATION*


if parm.neigh; k.sb.mask_ini = k.sb.mask; end

if parm.par
    delete(gcp('nocreate')); poolobj=parpool(parm.par_n);
    par_n_realisation = ceil(parm.n_realisation/poolobj.NumWorkers);
    
    YY=cell(poolobj.NumWorkers,1);
    tt=cell(poolobj.NumWorkers,1);
    
    parm_pool=parm;
    parm_pool.n_realisation = par_n_realisation;
    
    parfor pool_i=1:poolobj.NumWorkers
        [YY{pool_i}, tt{pool_i}]=BSGS_par_in(X, Z, kernel, k, Nscore, grid, parm_pool);
    end
    delete(poolobj)
   
   Y=YY{1};
   for pool_i=2:numel(YY)
       for i_scale=1:parm.n_scale
           Y{i_scale}.m = [Y{i_scale}.m; YY{pool_i}{i_scale}.m];
           Y{i_scale}.m_ns = [Y{i_scale}.m_ns; YY{pool_i}{i_scale}.m_ns];
       end
   end
    
    
else
    [Y,t.scale]=BSGS_par_in(X, Z, kernel, k, Nscore, grid, parm);
end


% save intial value
if parm.neigh; k.sb.mask = k.sb.mask_ini; end
clear X_ini k.sb.mask_ini


t.global = toc(t.global );


%% * 5. *SAVE IT*
if parm.saveit
    mkdir(['result/', parm.familyname])
    filename=['result/', parm.familyname, 'SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
    save(filename, 'parm', 'Y', 'grid', 't', 'X', 'Z', 'X_true', 'k', 'kernel', 'Nscore')
end


end