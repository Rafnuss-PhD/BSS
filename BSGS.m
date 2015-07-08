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

function [Y, U]=BSGS(X,Z,~,grid,plotit,n_sim)
% force input in column
X.x=X.x(:);X.y=X.y(:);Z.x=Z.x(:);Z.y=Z.y(:);

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
% * *KRIGING INPUT*

model     = 4; % model type
k.wradius = 1.3;
k.rotation  = 0;
k.range     = [160; 15]; % measure in unit

k.model = [model k.range(1) k.range(2) k.rotation; 1 1 1 1];
k.var   = [.99; 0.01];

k.nb_neigh  = [4 4 4 4; 15 15 15 15]; % min and max number of point


%%
% * *SEARCHING WINDOWS:* Super Block Grid
% creationg of the superblock grid
nx = 20; % number of superblock grid
ny = 20;
nb_max =50; % max number of point to take in the mask (random sampling)
plotit = 0; % plot
[k, X] = SuperBlockGridCreation(k, nx, ny, grid{end}.x(end), grid{end}.y(end), X, nb_max, plotit);
clear nx ny nb_max plotit


%% Non-parametric Relationship
% Link primary and secondary data.
kernel = kernel_est(X,Z);


%%
% * *NSCORE*
Nscore = nscore_perso(X.d,'linear',kernel);


%%
% * *SECONDARY*
scale=numel(grid);
Y=cell(scale,1);
U=cell(scale);
dy=cell(scale,1);

plotit=0;
if plotit
    figure('units','normalized','outerposition',[0 0 1 1]);
end

for s=1:scale % for each scale
    % Allocating space for resulting field. The third dimension is for at each simulation because the path is randomized and therefore the order of Y.x,Y.y and Y.d change.
    Y{s}.x=grid{s}.x; Y{s}.y=grid{s}.y; Y{s}.X=grid{s}.X; Y{s}.Y=grid{s}.Y; Y{s}.nx=grid{s}.nx; Y{s}.ny=grid{s}.ny;
    Y{s}.m=repmat({nan(grid{s}.ny,grid{s}.nx)},n_sim,1); % matrix des resutlats
    
    % populate the grid from previous scale.
    if s~=1
        for i_sim=1:n_sim
            Y{s}.m{i_sim}( 1:(grid{s-1}.dy/grid{s}.dy):end, 1:(grid{s-1}.dx/grid{s}.dx):end) = Y{s-1}.m{i_sim};
        end
    end
    
    % Assimilate the hard data (X) into the grid
    hard_data_idx=find(ismember(X.y,grid{s}.Y)&ismember(X.x,grid{s}.X));
    for i_sim=1:n_sim
        for hd=1:numel(hard_data_idx)
            Y{s}.m{i_sim}(X.x(hard_data_idx(hd))==grid{s}.X & X.y(hard_data_idx(hd))==grid{s}.Y) = X.d(hard_data_idx(hd));
        end
    end
    % remove the assimilated data from X.
    X.d(hard_data_idx)=[];
    X.x(hard_data_idx)=[];
    X.y(hard_data_idx)=[];
    k.sb.mask(:,:,hard_data_idx)=[];
    X.n=numel(X.d);
    clear hd i_sim hard_data_idx
    
    
    %% Generate the order for visting cells
    % Randomly permute the cell not known (to be visited). And generate the
    % Y.x and Y.y coordinate in a random order.
    % ATTENTION: the order is the same for all simulation but could be different
    rng('shuffle')
    Y{s}.sim.xy=grid{s}.xy(isnan(Y{s}.m{1}));
    Y{s}.sim.n=numel(Y{s}.sim.xy);
    Y{s}.sim.xy_r=Y{s}.sim.xy(randperm(Y{s}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
    [Y{s}.sim.x_r,Y{s}.sim.y_r] = ind2sub([grid{s}.ny, grid{s}.nx],Y{s}.sim.xy_r); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
    
    %%
    % * *RANDOM NORMAL FIELD* This create the random Normal distribution used for sampling the posteriori distribution at each point.
    U{s} = randn(Y{s}.sim.n,n_sim); % random field
    
    
    %% Simulation Loop
    % Here start stopeach simulation loop
    [Y{s}.m, dy{s}] = BSGS_pt(X, Z, Y{s}, normcdf(U{s}), kernel, Nscore, k, grid{s});
    
    if plotit
        subplot(scale,1,s);  hold on;axis tight
        imagesc(Y{s}.x,Y{s}.y,Y{s}.m{end})
        %scatter(Y{s}.X(:),Y{s}.Y(:),Y{s}.m(:))
        scatter(X.x,X.y,[],X.d,'o','filled','MarkerEdgeColor','k')
    end
    
    t(s+1)=toc;
    disp(['Simulation finish in : ' num2str(t(s+1)-t(s))])
end


end