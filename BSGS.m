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

function [Res, t, kern, k, parm, Nscore, filename] = BSGS(Prim,Sec,grid_gen,parm)
t.global = tic;
addpath('./kriging;./kernel;')

%% * *INPUT CEHCKING*
% This section of the code generates a valid parm structure based on the
% inputted parm. If some parm value haven't been input, this section will
% fill automatically with defautl value. This may not allowed be the best.

% Force input in column
Prim.x=Prim.x(:);Prim.y=Prim.y(:);Prim.d=Prim.d(:);
Sec.x=Sec.x(:);Sec.y=Sec.y(:);

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

% Path
if ~isfield(parm, 'path'),          parm.path            = 'linear'; end
if ~isfield(parm, 'path_random'),   parm.path_random     = true; end

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

% Kriging parameter
parm.k.covar = kriginginitiaite(parm.k.covar);
if ~isfield(parm, 'nscore'),        parm.nscore        =1; end % use normal score (strongly advice to use it.)
if ~isfield(parm, 'k') || ~isfield(parm.k, 'method'),  parm.k.method = 'smart'; end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'quad'),  parm.k.quad = 0; end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'nb'),  parm.k.nb = [0 0 0 0 0; 5 5 5 5 5]; end

if strcmp(parm.k.method,'sbss') && ~isfield(parm, 'k') || ~isfield(parm.k, 'sb') || ~isfield(parm.k.sb, 'nx') || ~isfield(parm.k.sb, 'ny') % super-block grid size (hard data)
    parm.k.sb.nx    = ceil(grid_gen.x(end)/parm.k.covar(1).range(1)*3);
    parm.k.sb.ny    = ceil(grid_gen.y(end)/parm.k.covar(1).range(2)*3);
end
if ~isfield(parm, 'k') || ~isfield(parm.k, 'wradius')
    parm.k.wradius  = Inf;
end
k = parm.k;

% Plot
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'bsgs'),  parm.plot.bsgs   =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'ns'),    parm.plot.ns     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'sb'),    parm.plot.sb     =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'kernel'),parm.plot.kernel =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'fitvar'),parm.plot.fitvar =0; end
if ~isfield(parm, 'plot') || ~isfield(parm.plot, 'krig'),  parm.plot.krig   =0; end

% Check the input for correct size dans dimension
assert(ismatrix(Sec.d),'Z is not 2Y.mat_simD');
assert(all([numel(Sec.y), numel(Sec.x)]==size(Sec.d)),'Z.x,y does not seems to be ok');
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
    
    grid{i_scale}.x=linspace(0, grid_gen.x(end), grid{i_scale}.nx)'; % coordinate of cells center
    grid{i_scale}.y=linspace(0, grid_gen.y(end), grid{i_scale}.ny)';
    grid{i_scale}.xy=(1:grid{i_scale}.nxy)';
    
    [grid{i_scale}.X, grid{i_scale}.Y] = meshgrid(grid{i_scale}.x,grid{i_scale}.y); % matrix coordinate
end


%% * 1. *SUPERBLOCK GRID CREATION*
% A mask (Boolean value) of the hard data is assigned to each superblock
% as follow: Only the n-closest (normalized by the covariance range) points
% (inside the ellipse/windows) to the centre of the superblock will be
% true. During the kriging, the mask of the superblock of the estimated
% point will be used to select the hard to add to the kriging system
if strcmp(parm.k.method,'sbss')
    k.sb.nx = parm.k.sb.nx; % number of superblock grid
    k.sb.ny = parm.k.sb.ny;
    [k, Prim] = SuperBlockGridCreation(k, grid_gen.x(end), grid_gen.y(end), Prim, parm.plot.sb, parm.k.nb(2,5));
end

%% * 2. *NON-PARAMETRIC RELATIONSHIP*
kern = kernel(Prim, Sec, parm.kernel,parm.plot.kernel);



%% * 3. *NORMAL SCORE TRANSFORM*
% Based on the hard data (well samples), a normal score transform is
% created using interpolation with power or exponential tail extrapolation.
% Matlab symbolique function are used for efficient coding. The back
% transform of the prior normal distribution function is also created
% (return the pdf in the initial space from the mean and variance in the
% normal space)
Nscore = nscore(kern, parm, parm.plot.ns); %Prim.d, kern.axis_prim, 'pchip', 'pchip', parm.plot.ns

% Create the normal space primary variable of known data
Prim.d_ns = Nscore.forward(Prim.d);





%% * 4. *RUN SIMULATION*

if strcmp(parm.k.method,'sbss'); k.sb.mask_ini = k.sb.mask; end % mask will be change after, we preserve its structure


if parm.par && parm.n_realisation > 1 % if parralelelisation is selected
    poolobj=gcp('nocreate');
    if isempty(poolobj)
        poolobj=parpool(parm.par_n);
    elseif poolobj.NumWorkers ~= parm.par_n
        delete(poolobj);
        poolobj=parpool(parm.par_n);
    end
    par_n_realisation = ceil(parm.n_realisation/poolobj.NumWorkers);
    
    RR=cell(poolobj.NumWorkers,1);
    tt=cell(poolobj.NumWorkers,1);

    parfor i_pool=1:poolobj.NumWorkers
        parm_pool=parm;
        parm_pool.n_realisation = par_n_realisation;
        parm_pool.i_pool = i_pool 
        [RR{i_pool}, tt{i_pool}]=BSGS_in(Prim, Sec, kern, k, Nscore, grid, parm_pool);
    end
    % delete(poolobj)
    
    Res=RR{1};    t.scale=[];    t.cstk=[];    t.pt=[];    t.krig=[];
    for i_pool=2:numel(RR)
        Res.m = [Res.m; RR{i_pool}.m];
        Res.m_ns = [Res.m_ns; RR{i_pool}.m_ns];
        Res.sim = [Res.sim; RR{i_pool}.sim];
        t.scale = [t.scale; tt{i_pool}.scale];
        t.cstk = [t.cstk; tt{i_pool}.cstk];
        t.pt = [t.pt; tt{i_pool}.pt];
        t.krig = [t.krig; tt{i_pool}.krig];
    end
    
    
elseif parm.n_realisation > 0
    [Res,t_out] = BSGS_in(Prim, Sec, kern, k, Nscore, grid, parm);
    t.scale = t_out.scale;
    t.cstk = t_out.cstk;
    t.pt = t_out.pt;
    t.krig = t_out.krig;
else
    Res={};
end



% save intial value
if strcmp(parm.k.method,'sbss'); k.sb.mask = k.sb.mask_ini; end



t.global = toc(t.global);


%% * 5. *SAVE IT*
filename=['result/', parm.familyname, 'SIM-', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'];
if parm.saveit
    mkdir(['result/', parm.familyname])
    save(filename, 'parm', 'Res', 'grid', 't', 'Prim', 'Sec', 'k', 'kern', 'Nscore')
end

if parm.notify
    unix(['echo "Simulation ' parm.familyname ' has finish now (' datestr(now,'yyyy-mm-dd_HH-MM-SS') ') in ' num2str(t.global/60) 'min" | mail -s "Simulation Finish" ' parm.notify_email]);
    % load handel; sound(y,Fs)
end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Res_out,t_out]=BSGS_in(Prim, Sec, kern, k, Nscore, grid, parm)

clear parm.k

Res = cell(parm.n_scale,1);
Res{end}.w=cell(parm.n_realisation,1);
t.toc.scale = [];
t.toc.cstk = [];
t.toc.pt = [];
t.toc.krig = [];


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
    
    % if strcmp(parm.k.method,'sbss'); k.sb.mask(:,:,hard_data_idx)=[]; end
    
    % Create the normal space result matrix and populate with known value
    Res{i_scale}.m_ns = repmat({nan(size(Res{i_scale}.m{1}))},parm.n_realisation,1);
    I=~isnan(Res{i_scale}.m{1});
    for i_realisation=1:parm.n_realisation
        Res{i_scale}.m_ns{i_realisation}(I) = Nscore.forward(Res{i_scale}.m{i_realisation}(I));
    end
    
    
    %% * 1.2. *SPIRAL SEARCH*
    % Create the windows for kringing with the function to compute the
    % normalized distence and the order of visit of the cells. Spiral
    % Search setting: previously data (on grid{i_scale} location)
    
    if strcmp(parm.k.method,'sbss')
        k.ss.el.dw = ceil(min(max(k.covar(1).range*k.wradius./[grid{i_scale}.dx grid{i_scale}.dy]), max(grid{i_scale}.nx,grid{i_scale}.ny)));
        
        [k.ss.el.X, k.ss.el.Y] = meshgrid(-k.ss.el.dw:k.ss.el.dw);% grid{i_scale} of searching windows
        [k.ss.el.X_T, k.ss.el.Y_T]=rotredtrans(k.ss.el.X*grid{i_scale}.dx, k.ss.el.Y*grid{i_scale}.dy, k.covar(1).azimuth, k.covar(1).range); % transforms the grid{i_scale}
        k.ss.el.dist = sqrt(k.ss.el.X_T.^2 + k.ss.el.Y_T.^2); % find distence
        [k.ss.el.dist_s, k.ss.el.dist_idx] = sort(k.ss.el.dist(:)); % sort according distence.
        k.ss.el.X_s=k.ss.el.X(k.ss.el.dist_idx); % sort the axis
        k.ss.el.Y_s=k.ss.el.Y(k.ss.el.dist_idx);
        k.ss.el.n = numel(k.ss.el.dist_idx);
    end
    
    %% * 1.3. *Generate the order for visting cells*
    % Randomly permute the cell not known (to be visited). And generate the
    % Y.x and Y.y coordinate in a random order.
    
    if i_scale<parm.cstk_s
        parm.cstk=0;
    else
        parm.cstk=1;
    end
    
    [Res{i_scale}.sim] = definepath(Res{i_scale},grid{i_scale},parm);
    
    
    %% * 1.4. *RANDOM NORMAL FIELD*
    % This create the random Normal distribution used for sampling the
    % posteriori distribution at each point.
    Res{end}.U = normcdf(randn(Res{i_scale}.sim.n,parm.n_realisation)); % random field
    
    
    %% * 2. *Simulation of Point*
    i_plot=0; % used for computing the number of point simulated. Used for ploting
    for i_pt=1:Res{i_scale}.sim.n; % loop over each point
        i_plot=i_plot+1;
        pt.i=i_pt;
        t.tic.pt = tic;
        
        if parm.cstk % if option constant weight is activate.
            % Find the current point position on the grid{i_scale}. It
            % will be used on each realisation.
            pt.x = Res{i_scale}.sim.x_r{1}(i_pt);
            pt.y = Res{i_scale}.sim.y_r{1}(i_pt);
            
            % Kriging system
            t.tic.krig = tic;
            pt.krig = kriging(pt,Res{i_scale},Prim,k,parm,1);
            t.toc.krig = [t.toc.krig; toc(t.tic.krig)];
            
            % Secondary Info
            % We first find the secondary value (and error) to create a
            % pdf. This pdf is then multiply inside the kernel to get the
            % density.
            % 1. Option using Sec.std
            % pt.sec.sec_pdf = normpdf(kern.axis_sec, Sec.d(Sec.y==Res{i_scale}.Y(pt.y,pt.x), Sec.x==Res{i_scale}.X(pt.y,pt.x))  , Sec.std(Sec.y==Res{i_scale}.Y(pt.y,pt.x), Sec.x==Res{i_scale}.X(pt.y,pt.x)));
            
            % 2. Option 
            pt.sec.sec_pdf = zeros(numel(kern.axis_sec),1);
            [~,id] = min( (kern.axis_sec-Sec.d(Sec.y==Res{i_scale}.Y(pt.y,pt.x), Sec.x==Res{i_scale}.X(pt.y,pt.x))).^2 );
            pt.sec.sec_pdf(id) = 1;
            
            pt.sec.dens = bsxfun(@times, kern.dens, pt.sec.sec_pdf');
            pt.sec.pdf = sum(pt.sec.dens, 2)/sum(pt.sec.dens(:));
        end
        
        
        %% * 3. *Simulation of Realisation*
        for i_realisation=1:parm.n_realisation
            if ~parm.cstk
                % Find the current point position on the grid{i_scale}.
                % This is changing for each realisation
                pt.x = Res{i_scale}.sim.x_r{i_realisation}(i_pt);
                pt.y = Res{i_scale}.sim.y_r{i_realisation}(i_pt);
                
                % Kriging
                t.tic.krig = tic;
                pt.krig = kriging(pt,Res{i_scale},Prim,k,parm,i_realisation);
                t.toc.krig = [t.toc.krig; toc(t.tic.krig)];
                
                % Secondary Info
                % 1. Option using Sec.std
                % pt.sec.sec_pdf = normpdf(kern.axis_sec, Sec.d(Sec.y==Res{i_scale}.Y(pt.y,pt.x), Sec.x==Res{i_scale}.X(pt.y,pt.x))  , Sec.std(Sec.y==Res{i_scale}.Y(pt.y,pt.x), Sec.x==Res{i_scale}.X(pt.y,pt.x)));

                % 2. Option 
                pt.sec.sec_pdf = zeros(kern.n,1);
                [~,id] = min( (kern.axis_sec-Sec.d(Sec.y==Res{i_scale}.Y(pt.y,pt.x), Sec.x==Res{i_scale}.X(pt.y,pt.x))).^2 );
                pt.sec.sec_pdf(id) = 1;

                pt.sec.dens = bsxfun(@times, kern.dens, pt.sec.sec_pdf);
                pt.sec.pdf = sum(pt.sec.dens, 2)/sum(pt.sec.dens(:));
            end
            
            if ~isempty(pt.krig.mask.prim) || ~isempty(pt.krig.mask.res)
                pt.krig.m = pt.krig.lambda'* [Prim.d_ns(pt.krig.mask.prim) ; Res{i_scale}.m_ns{i_realisation}(pt.krig.mask.res)];
            else % no neighbor
                pt.krig.m=0; % mean ?
            end
            
            %% * 3.1 *KRIGING DISTRIBUTION*
            % Back transform the normal distribution (from krigeage) in
            % original space in the  grid{i_scale}. this become
            % the prior distribution
            pt.krig.pdf = Nscore.dist(pt.krig.m, sqrt(pt.krig.s));
            pt.krig.pdf = pt.krig.pdf./sum(pt.krig.pdf);
            
            
            %% * 3.2 *SECONDARY DISTRIBUTION*
            % Multiply the (nomalized) prior and likelihood to get the
            % (normalised) posteriori
            %pt.sec.pdf(pt.sec.pdf==0)=eps; % if sec and krig don't have overlaps, put eps.
            %pt.krig.pdf(pt.krig.pdf==0)=eps;
            % pt.aggr.pdf = kern.prior.^(1-parm.w_krig(i_scale)-parm.w_sec(i_scale)) .* pt.sec.pdf.^parm.w_sec(i_scale) .* pt.krig.pdf.^parm.w_krig(i_scale);
            w = aggr_fx(Res,Sec,parm,grid,i_realisation,i_scale,pt);
            Res{end}.w{i_realisation}=[Res{end}.w{i_realisation};w];

            % pt.aggr.pdf = pt.sec.pdf.^(1-w) .* pt.krig.pdf.^w;
            % pt.aggr.pdf = pt.sec.pdf.^1 .* pt.krig.pdf.^1;
            % pt.aggr.pdf = kern.prior.^(1-parm.aggr.sum).* pt.sec.pdf.^(parm.aggr.sum-w) .* pt.krig.pdf.^w; 
            pt.aggr.pdf = kern.prior.^(w-1).* pt.sec.pdf.^(1-w) .* pt.krig.pdf.^1; 
            pt.aggr.pdf(isnan(pt.aggr.pdf) | pt.aggr.pdf==Inf)=0;


            pt.aggr.pdf = pt.aggr.pdf./sum(pt.aggr.pdf);
            
            
            %% * 3.3  *SAMPLING*:
            % Sample a point in the posteri distribution. To do so the CDF
            % is created and interpolation tool is used to assure
            % randomized sampling.
            pt.aggr.cdf = cumsum(pt.aggr.pdf);
            
            % Interpolation only works if it is mono.. increasing. But the
            % cdf can reach 1 so we just add eps (very small value).
            if ~all(diff(pt.aggr.cdf)>0)
                pt.aggr.cdf = pt.aggr.cdf +  linspace(0,numel(pt.aggr.cdf)*eps*2,numel(pt.aggr.cdf))';
                assert(all(diff(pt.aggr.cdf)>0))
            end
            
            pt.sampled = interp1(pt.aggr.cdf, kern.axis_prim, Res{end}.U(i_pt,i_realisation),'pchip');
            
            
            %% * 3.4 *PLOTIT*
            if parm.plot.krig && i_realisation==1 && (i_plot==1 || mod(i_plot+49,50)==0  || Nscore.forward(pt.sampled)<-3  || Nscore.forward(pt.sampled)>3 ) %
                figure(1); clf
                
                subplot(3,2,[1 4]);
                hold on
                h1=imagesc(Res{i_scale}.x,Res{i_scale}.y,Res{i_scale}.m_ns{1},'AlphaData',~isnan(Res{i_scale}.m_ns{1}));
                
                if strcmp(parm.k.method,'sbss');
                    sb_i = min([round((Res{i_scale}.y(pt.y)-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny]);
                    sb_j = min([round((Res{i_scale}.x(pt.x) -k.sb.x(1))/k.sb.dx +1)'; k.sb.nx]);
                    windows=nan(k.sb.ny,k.sb.nx);
                    for u=1:length(k.el_X) %.. look at all point...
                        if sb_i+k.el_X(u)<=k.sb.nx && sb_j+k.el_Y(u)<=k.sb.ny && sb_i+k.el_X(u)>=1 && sb_j+k.el_Y(u)>=1% check to be inside the grid
                            windows(sb_i+k.el_X(u), sb_j+k.el_Y(u))=1;
                        end
                    end
                    h2=imagesc(k.sb.x,k.sb.y,windows,'AlphaData',windows*.5);
                    h3=mesh([0 k.sb.x+k.sb.dx/2],[0 k.sb.y+k.sb.dy/2],zeros(k.sb.ny+1, k.sb.nx+1),'EdgeColor','k','facecolor','none');
                end
                
                tt=-pi:0.01:pi;
                x=Res{i_scale}.x(pt.x)+k.covar.range(1)*cos(tt);
                x2=Res{i_scale}.x(pt.x)+k.wradius*k.covar.range(1)*cos(tt);
                y=Res{i_scale}.y(pt.y)+k.covar.range(2)*sin(tt);
                y2=Res{i_scale}.y(pt.y)+k.wradius*k.covar.range(2)*sin(tt);
                h4=plot(x,y,'--r'); h5=plot(x2,y2,'-r');
                h6=plot([Res{i_scale}.x(pt.x) Res{i_scale}.x(pt.x)],[min(y2) max(y2)],'-r');
                h7=plot([min(x2) max(x2)], [Res{i_scale}.y(pt.y) Res{i_scale}.y(pt.y)],'-r');
                
                lambda_c= 36+60.*(abs(pt.krig.lambda)-min(abs(pt.krig.lambda)))./range(abs(pt.krig.lambda));
                
                h8=scatter(Prim.x,Prim.y,[],Prim.d_ns,'s','filled');
                
                n_hd = numel(Prim.x(pt.krig.mask.prim));
                sel=[Prim.x(pt.krig.mask.prim) Prim.y(pt.krig.mask.prim); Res{i_scale}.X(pt.krig.mask.res) Res{i_scale}.Y(pt.krig.mask.res)];
                XY_ns = [Prim.d_ns(pt.krig.mask.prim) ; Res{i_scale}.m_ns{i_realisation}(pt.krig.mask.res)];
                h9=scatter(sel(1:n_hd,1),sel(1:n_hd,2),lambda_c(1:n_hd),XY_ns(1:n_hd),'s','filled','MarkerEdgeColor','k');
                h10=scatter(sel(n_hd+1:end,1),sel(n_hd+1:end,2),lambda_c(n_hd+1:end),XY_ns(n_hd+1:end),'o','filled','MarkerEdgeColor','k');
                h11=scatter(Res{i_scale}.x(pt.x), Res{i_scale}.y(pt.y),100,pt.krig.lambda'* XY_ns(:),'o','filled','MarkerEdgeColor','r','LineWidth',1.5);
                
                xlabel('x[m]');ylabel('y[m]');
                %colorbar;
                xlim([Res{i_scale}.x(1) Res{i_scale}.x(end)])
                ylim([Res{i_scale}.y(1) Res{i_scale}.y(end)])
                %set(gca,'YDir','reverse');
                
                % legend([h3 h6 h8 h9 h10 h11],'Super grid','Window search with quadrant','Hard data point','Selected hard data point','Selected previously simulated point','Simulated Point','Location','northoutside','Orientation','horizontal')
                
                
                subplot(3,2,5); hold on;
                plot( kern.axis_prim,pt.krig.pdf)
                plot( kern.axis_prim,pt.sec.pdf)
                plot( kern.axis_prim,pt.aggr.pdf)
                plot(  kern.axis_prim, max(pt.aggr.pdf)*pt.aggr.cdf)
                plot([pt.sampled pt.sampled],[0 max(pt.aggr.pdf)],'k')
                legend(['Kriging, w=' num2str(w)], ['Secondary, w=' num2str(1-w)],'Aggregated','Aggreag. cdf','sampled location')
                
                
                subplot(3,2,6); hold on;
                %ksdensity(Sec.d(:));
                %ksdensity(Res{i_scale}.m{i_realisation}(~isnan(Res{i_scale}.m{i_realisation}(:))));
                %[f,x]=hist(X_ini.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
                [f,x]=hist(Prim.d_ns(:)); plot(x,f/trapz(x,f),'linewidth',2);
                [f,x]=hist(Res{i_scale}.m_ns{i_realisation}(:),50); plot(x,f/trapz(x,f));
                legend('Well sampling', 'Simulation')
                
                drawnow
                keyboard
            end
            
            %% 3.5 *Back transform*
            Res{i_scale}.m{i_realisation}(pt.y,pt.x) = pt.sampled;
            Res{i_scale}.m_ns{i_realisation}(pt.y,pt.x) = Nscore.forward(pt.sampled); % add only the last simulated point to normal score
        end
        % NOTHING SHOULD BE DONE AFTER THAT WHICH INVOLVED THE WEIGHT
        % BECAUSE WE CHANGED RES.M...
    end

end
t_out = t.toc;

Res_out = Res{end};
Res_out.sim = [];
%Res_out.U = [];
for i_scale=1:parm.n_scale
    Res_out.sim = [Res_out.sim Res{i_scale}.sim];
    %Res_out.sim = [Res_out.U Res{i_scale}.U];
end



end


function w = aggr_fx(Res,Sec,parm,grid,i_realisation,i_scale,pt)


i_t = mod(i_realisation,size(parm.aggr.T,1));
if i_t==0;
    i_t=size(parm.aggr.T,1);
end
x = pt.i./ (grid{end}.nxy - Res{i_scale}.nxy+Res{i_scale}.sim.n);
assert(x<=1,'error')

switch parm.aggr.method
    case 'cst'
        w=parm.aggr.T(i_t);
    case 'step'
        if (x<parm.aggr.T(i_t))
            w=0;
        else 
            w=1;
        end
    case 'linear'
        if (x<parm.aggr.T(i_t,1))
            w  =0;
        elseif (x>parm.aggr.T(i_t,2))
            w = 1;
        else 
            w = ( parm.aggr.T(i_t,2)-x ) / (parm.aggr.T(i_t,2)-parm.aggr.T(i_t,1));
        end
    case 'sigmoid'
        a = parm.aggr.T(i_t,1);
        b = parm.aggr.T(i_t,2);
        w = (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
    otherwise
        error('no Aggr method defined')  
end

% w = parm.aggr.w(i_realisation);

end
