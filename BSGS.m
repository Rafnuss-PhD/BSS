%% BSS is computing the Bayesian Sequential Gaussian Simulation.
% It's a sequential stochatic simulation of a primary variable from
% combination of sparce primary (X) data and more dence secondary data (Z).
% The grid is define by the size of the secondary variable (Z). The primary
% variable can be at any point of this grid (X_x,X_y)
%
% INPUT:
%
% * X.d     : Primary variable (1D-X.n elts)
% * X.x,y   : x,y-coordinate of the primary variable (1D-X.n elts)
% * X.y     : y-coordinate of the primary variable (1D-X.n elts)
% * Z.d     : Secondary variable (2D- grid.nx x grid.n elts)
% * Z.x,y   : x,y-coordinate of the secondary variable (1D elts)
% * Z.std   : Secondary variable std error (2D- grid.nx x grid.n elts)
%
% OUTPUT:
%
% * Y_mat   : Simulated Primary variable in matrix (3rd dim for simulation)
% * U       :
% * t_sim   : time of simulation
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

function [Y, U]=BSGS(X,Z,Y_true,grid)
% force input in column
X.x=X.x(:);X.y=X.y(:);Z.x=Z.x(:);Z.y=Z.y(:);
figure;
%% Input checkout
% Check the input for correct size dans dimension
assert(ismatrix(Z.d),'Z is not 2Y.mat_simD');
assert(all([numel(Z.y), numel(Z.x)]==size(Z.d)),'Z.x,y does not seems to be ok');
assert(size(X.d,2)==1,'X.d is not a vertical vector 1D');
assert(size(X.x,2)==1,'X.dx is not a vertical vector 1D');
assert(size(X.y,2)==1,'X.y is not a vertical vector 1D');
assert(all(size(X.y)==size(X.x)),'X.x and X.y don''t have the same dimension');
assert(all(size(X.d)==size(X.x)),'X.d and X.x (or X.y) don''t have the same dimension');

%% Allocate space and setting basic information.

%%
% * *KRIGING input*

% [xrange, yrange, C0]=fit_variogramm(Y,grid{end});

k.model     = 'exp'; % model type
k.var       = 1; % data are normalized later on...
k.range     = [160 15]; % measure in unit
k.wradius   = 1.4; % as a factor of x and yrange
k.nb_neigh  = [4 4 4 4; 15 15 15 15]; % min and max number of point
k.ang       = 0;
k.C_0       = 0.1; 

switch k.model % portee or radius is remove because dist in normalized to x,y range.
    case 1
    case 2
    case 'exp'
        k.C_fx = @(h) k.var*exp( -3* h ) + k.C_0;
    case 'sph'
        k.C_fx = @(h)  k.var*(1-1.5*h + .5*h.^3);
        %k.gamma_fx =@(h) k.var*(1.5*h - .5*h.^3);
end

%%
% * *SEARCHING WINDOWS:* Super Block Grid
% creationg of the superblock grid
nx = 20;
ny = 20;
[k, X] = SuperBlockGridCreation(k, nx, ny, grid{end}.x(end), grid{end}.y(end), X);
clear nx ny

%% Non-parametric Relationship
% Link primary and secondary data.
kernel = kernel_est(X,Z);


%%
% * *SECONDARY*
scale=numel(grid);
Y=cell(scale,1);
U=cell(scale,1);
dy=cell(scale,1);

for s=1:scale % for each scale
    % Alocating space for resulting field. The third dimension is fo at each simulation because the path is randomized and therefore the order of Y.x,Y.y and Y.d change.
    Y{s}.x=grid{s}.x; Y{s}.y=grid{s}.y; Y{s}.X=grid{s}.X; Y{s}.Y=grid{s}.Y;
    Y{s}.m=nan(grid{s}.ny,grid{s}.nx); % matrix des resutlats
    
    % populate the grid from previous scale.
    if s~=1
        Y{s}.m( 1:(grid{s-1}.dy/grid{s}.dy):end, 1:(grid{s-1}.dx/grid{s}.dx):end) = Y{s-1}.m;
    end
    
    % Assimilate the hard data (X) into the grid
    hard_data_idx=find(ismember(X.y,grid{s}.Y)&ismember(X.x,grid{s}.X));
    for hd=1:numel(hard_data_idx)
        Y{s}.m(X.x(hard_data_idx(hd))==grid{s}.X & X.y(hard_data_idx(hd))==grid{s}.Y)=X.d(hard_data_idx(hd));
    end
    % remove the assimilated data from X.
    X.d(hard_data_idx)=[];
    X.x(hard_data_idx)=[];
    X.y(hard_data_idx)=[];
    k.sb.mask(:,:,hard_data_idx)=[];
    X.n=numel(X.d);
    clear hd
    
    
    %% Generate the order for visting cells
    % Randomly permute the cell not known (to be visited). And generate the
    % Y.x and Y.y coordinate in a random order.
    % ATTENTION: the order is the same for all simulation but could be different
    rng('shuffle')
    Y{s}.sim.xy=grid{s}.xy(isnan(Y{s}.m));
    Y{s}.sim.n=numel(Y{s}.sim.xy);
    Y{s}.sim.xy_r=Y{s}.sim.xy(randperm(Y{s}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
    [Y{s}.sim.x_r,Y{s}.sim.y_r] = ind2sub([grid{s}.ny, grid{s}.nx],Y{s}.sim.xy_r); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
    
    %%
    % * *RANDOM NORMAL FIELD* This create the random Normal distribution used for sampling the posteriori distribution at each point.
    U{s} = randn(Y{s}.sim.n,1); % random field
    
    
    %% Simulation Loop
    % Here start stopeach simulation loop
    [Y{s}.m, dy{s}] = BSGS_pt(X, Z, Y{s}, normcdf(U{s}), kernel, k, grid{s});
    
    if s<scale
        subplot(scale-1,1,s);
    else
        figure;
    end
    hold on
    imagesc(Y{s}.x,Y{s}.y,Y{s}.m)
    %scatter(Y{s}.X(:),Y{s}.Y(:),Y{s}.m(:))
    scatter(X.x,X.y,[],X.d,'o')
    xlim([0 Y{s}.x(end)])
    ylim([0 Y{s}.y(end)])
    axis equal
    
    t(s+1)=toc;
    disp(['Simulation finish in : ' num2str(t(s+1)-t(s))])
end

a=round(sqrt(scale));
b=ceil(scale/a);
figure;
for i=1:scale
    subplot(a,b,i);
    hist(dy{i})
end

end