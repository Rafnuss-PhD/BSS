%% BSS is computing the Bayesian Sequential Gaussian Simulation.
% It's a sequential stochatic simulation of a primary variable from
% combination of sparce primary (X) data and more dence secondary data (Z).
% The grid is define by the size of the secondary variable (Z). The primary
% variable can be at any point of this grid (X.x,X.y)
%
% INPUT:
%
% * X.d     : Primary variable (1D-X.n elts)
% * X.x,y   : x,y-coordinate of the primary variable (1D-X.n elts)
% * Z.d     : Secondary variable (2D- grid.nx x grid.n elts)
% * Z.x,y   : x,y-coordinate of the secondary variable (1D elts)
% * Z.std   : Secondary variable std error (2D-grid.nx x grid.n elts)
%
% OUTPUT:
%
% * Y       : Simulated Primary variable in matrix
% * t       : Time of simulations
% * kernel  : kernel information
% * k     
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

function [Y, t, kernel, k] = BSGS(X,Z,X_true,grid,parm)
t.tic.global = tic;
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


%%
% * *VARIOGRAM INPUT*

if parm.fitvar
    [k.range, fig.varfit] = fit_variogramm(X,Z,parm.plot.fitvar);
    disp(['The fitted range is: ', num2str(k.range)])
    keyboard
else % default variogramm 
    k.range = [100 10];
end
k.rotation = 0;
k.model = [4 k.range k.rotation; 1 1 1 1];
k.wradius = 1.3;
k.var   = [.99; 0.01];
k.nb_neigh  = [5 5 5 5; 30 30 30 30]; % min and max number of point

%%
% * *SEARCHING WINDOWS:* 
% creationg of the Super Block Grid 
if parm.neigh
    nx = 20; % number of superblock grid
    ny = 20;
    nb_max = 40; % max number of point to take in the mask (random sampling)
    [k, X] = SuperBlockGridCreation(k, nx, ny, grid{end}.x(end), grid{end}.y(end), X, nb_max, parm.plot.sb);
    clear nx ny nb_max plotit2
end

%% Non-parametric Relationship
% Link primary and secondary data.
kernel = kernel_est(X, Z, parm.plot.kernel);

%%
% * *NSCORE*
Nscore = nscore_perso(X.d, 'linear', 'linear', kernel, parm.plot.ns);
kernel.y_ns = Nscore.forward(kernel.y); % convert the kernel grid in normal space
[~,kernel.y_ns_unique,~] = unique(kernel.y_ns); % find the unique value as during transformation some become very similar

%%
% * *Scale to simulate*
assert(max(parm.scale)<=numel(grid))
Y=cell(numel(parm.scale),1);
X.x_ini = X.x; X.y_ini = X.y; X.d_ini = X.d;

if parm.neigh; k.sb.mask_ini = k.sb.mask; end

for scale_i=1:numel(parm.scale) % for each scale
    t.tic.scale = tic;
    s = parm.scale(scale_i);
    % Allocating space for resulting field. The third dimension is for at each simulation because the path is randomized and therefore the order of Y.x,Y.y and Y.d change.
    Y{scale_i}.x=grid{s}.x; 
    Y{scale_i}.y=grid{s}.y; 
    Y{scale_i}.X=grid{s}.X; 
    Y{scale_i}.Y=grid{s}.Y; 
    Y{scale_i}.nx=grid{s}.nx; 
    Y{scale_i}.ny=grid{s}.ny;
    Y{scale_i}.m=repmat({nan(grid{s}.ny,grid{s}.nx)},parm.n_realisation,1); % matrix des resutlats
    
    % Populate the grid from previous scale.
    if scale_i~=1 % not at first scale
        for i_sim=1:parm.n_realisation
            Y{scale_i}.m{i_sim}( 1:(grid{s-1}.dy/grid{s}.dy):end, 1:(grid{s-1}.dx/grid{s}.dx):end) = Y{s-1}.m{i_sim};
        end
    end
    
    % Assimilate the hard data (X) into the grid
    hard_data_idx=find(ismember(X.y,grid{s}.Y)&ismember(X.x,grid{s}.X));
    for i_sim=1:parm.n_realisation
        for hd=1:numel(hard_data_idx)
            Y{scale_i}.m{i_sim}(X.x(hard_data_idx(hd))==grid{s}.X & X.y(hard_data_idx(hd))==grid{s}.Y) = X.d(hard_data_idx(hd));
        end
    end
    
    % Remove the assimilated data from X.
    X.d(hard_data_idx)=[];
    X.x(hard_data_idx)=[];
    X.y(hard_data_idx)=[];
    if parm.neigh; k.sb.mask(:,:,hard_data_idx)=[];  end
    X.n=numel(X.d);
    clear hd i_sim hard_data_idx
    
    % Simulation
    [Y{scale_i}.m, Y{scale_i}.m_ns] = BSGS_s(X, Z, Y{scale_i}, kernel, Nscore, k, grid{s}, parm);
   
    % display info
    disp(['Simulation ' num2str(scale_i) '/' num2str(numel(parm.scale)) ' finished in : ' num2str(toc(t.tic.scale))])
    t.scale{scale_i} = toc(t.tic.scale);
end

X.x = X.x_ini; X.y = X.y_ini; X.d = X.d_ini; X.n=numel(X.d);
if parm.neigh; k.sb.mask = k.sb.mask_ini; end
clear X.d_ini X.y_ini X.x_ini k.sb.mask_ini

t.global = toc(t.tic.global);
clear t.tic


%% Save it
if parm.saveit
    save(['result/', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM'), '.mat'], 'parm', 'Y', 'grid', 't', 'X', 'Z', 'X_true', 'k', 'kernel', 'Nscore')
end

end