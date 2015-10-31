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
% * parm    : Parameters of the simulation
%
% OUTPUT:
%
% * Y       : Simulated Primary variable
% * t       : Time of simulations
% * kernel  : kernel information
% * k       : kriging information
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

function [Y, t, kernel, k] = BSGS(X,Z,X_true,grid,parm)
t.tic.global = tic;

%% 
% * *INPUT CHECKOUT*
% force input in column
X.x=X.x(:);X.y=X.y(:);Z.x=Z.x(:);Z.y=Z.y(:);
% Check the input for correct size dans dimension
assert(ismatrix(Z.d),'Z is not 2Y.mat_simD');
assert(all([numel(Z.y), numel(Z.x)]==size(Z.d)),'Z.x,y does not seems to be ok');
assert(size(X.d,2)==1,'X.d is not a vertical vector 1D');
assert(size(X.x,2)==1,'X.dx is not a vertical vector 1D');
assert(size(X.y,2)==1,'X.y is not a vertical vector 1D');
assert(all(size(X.y)==size(X.x)),'X.x and X.y don''t have the same dimension');
assert(all(size(X.d)==size(X.x)),'X.d and X.x (or X.y) don''t have the same dimension');

%%
% * *DEFAULT VALUE*
if ~isfield(parm, 'seed'),   parm.seed       =rand(); end
if ~isfield(parm, 'covar'),  parm.covar      =parm.gen.covar; end
if ~isfield(parm, 'saveit'), parm.saveit     =1; end
if ~isfield(parm, 'name'),   parm.name       =''; end
if ~isfield(parm, 'unit'),   parm.unit       =''; end
if ~isfield(parm, 'n_realisation'),   parm.n_realisation       =1; end
% Run option
if ~isfield(parm, 'neigh'),  parm.neigh      =1; end
if ~isfield(parm, 'nscore'), parm.nscore     =1; end
if ~isfield(parm, 'cstk'),   parm.cstk       =1; end
if ~isfield(parm, 'fitvar'), parm.fitvar     =0; end
% Plot
if isfield(parm, 'plot')
    if ~isfield(parm.plot, 'bsgs'),  parm.plot.bsgs   =0; end
    if ~isfield(parm.plot, 'ns'),    parm.plot.ns     =0; end
    if ~isfield(parm.plot, 'sb'),    parm.plot.sb     =0; end
    if ~isfield(parm.plot, 'kernel'),parm.plot.kernel =0; end
    if ~isfield(parm.plot, 'fitvar'),parm.plot.fitvar =0; end
    if ~isfield(parm.plot, 'krig'),  parm.plot.krig   =0; end
else
    parm.plot.bsgs   =0;
    parm.plot.ns     =0;
    parm.plot.sb     =0;
    parm.plot.kernel =0;
    parm.plot.fitvar =0;
    parm.plot.krig   =0;
end
% Kriging parameter
if ~isfield(parm, 'nb_neigh'),   parm.nb_neigh    = [2 2 2 2 0; 7 7 7 7 10]; end
if isfield(parm, 'k')
    if isfield(parm.k, 'range')
        if ~isfield(parm.k.range, 'min'),parm.k.range.min = [min(X.d(:))-2 min(X.d(:))-2]; end
        if ~isfield(parm.k.range, 'max'),parm.k.range.max = [max(X.d(:))+2 max(X.d(:))+2]; end
    else
        parm.k.range.min = [min(X.d(:))-2 min(X.d(:))-2];
        parm.k.range.max = [max(X.d(:))+2 max(X.d(:))+2];
    end
    if isfield(parm.k, 'sb')
        if ~isfield(parm.k.sb, 'nx'),    parm.k.sb.nx     = ceil(grid{end}.x(end)/parm.covar.modele(1,2)*3); end
        if ~isfield(parm.k.sb, 'ny'),    parm.k.sb.ny     = ceil(grid{end}.y(end)/parm.covar.modele(1,3)*3); end
    else
        parm.k.sb.nx     = ceil(grid{end}.x(end)/parm.covar.modele(1,2)*3);
        parm.k.sb.ny     = ceil(grid{end}.y(end)/parm.covar.modele(1,3)*3);
    end
else
    parm.k.range.min = [min(X.d(:))-2 min(X.d(:))-2];
    parm.k.range.max = [max(X.d(:))+2 max(X.d(:))+2];
    parm.k.sb.nx     = ceil(grid{end}.x(end)/parm.covar.modele(1,2)*3);
    parm.k.sb.ny     = ceil(grid{end}.y(end)/parm.covar.modele(1,3)*3);
end


%%
% * *VARIOGRAM INPUT*
if parm.fitvar
    [k.range, k.var] = fit_variogramm(X,Z,parm.plot.fitvar, X_true);
else % default variogramm 
    k.range  = parm.covar.modele(1,2:3);
    k.var    = parm.covar.c;
end
k.rotation = 0;
k.model = [4 k.range k.rotation; 1 1 1 1];
k.wradius = 1.3;
k.nb_neigh  = parm.nb_neigh; % min and max number of point


%%
% * *SEARCHING WINDOWS:* 
% creationg of the Super Block Grid 
if parm.neigh
    k.sb.nx = parm.k.sb.nx; % number of superblock grid
    k.sb.ny = parm.k.sb.ny;
    [k, X] = SuperBlockGridCreation(k, grid{end}.x(end), grid{end}.y(end), X, parm.plot.sb,parm.nb_neigh(2,5));
end

%% Non-parametric Relationship
% Link primary and secondary data.
kernel = kernel_est(X, Z, parm.k.range, parm.plot.kernel);


%%
% * *NSCORE*
if parm.nscore
    Nscore = nscore_perso(X.d, 'linear', 'linear', kernel, parm.plot.ns);
else
    Nscore.forward = @(x) x;
    Nscore.inverse = @(x) x;
    Nscore.dist    = @(mu,sigma) normpdf(kernel.y,mu,sigma)/sum(normpdf(kernel.y,mu,sigma)); 
end

%%
% * *INITIATE*
assert(max(parm.scale)<=numel(grid)) % assert that simulation scale exist
Y=cell(numel(parm.scale),1); % allocate variable space
X_ini = X; %keep initial data intact

if parm.neigh; 
    k.sb.mask_ini = k.sb.mask; 
end

for scale_i=1:numel(parm.scale) % for each scale
    t.tic.scale = tic;
    s = parm.scale(scale_i); % find scale level
    
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
    if parm.neigh; 
        k.sb.mask(:,:,hard_data_idx)=[];
    end
    X.n=numel(X.d);
    clear hd i_sim hard_data_idx
    
    %%
    % * *RESULT ALLOCATION*
    % Create the normal space result matrix and populate with known value 
    Y{scale_i}.m_ns = repmat({nan(size(Y{scale_i}.m{1}))},parm.n_realisation,1); 
    I=~isnan(Y{scale_i}.m{1});
    for i_realisation=1:parm.n_realisation
        Y{scale_i}.m_ns{i_realisation}(I) = Nscore.forward(Y{scale_i}.m{i_realisation}(I));
    end

    % Create the normal space primary variable of known data
    X.d_ns = Nscore.forward(X.d);
    
    %%
    % * *SEARCHING-WINDOWS*: create the windows for kringing with the function to compute
    % the normalized distence and the order of visit of the cells.
    % Spiral Search setting: previously data (on grid{s} location)
    if parm.neigh
        [k.ss.el.X, k.ss.el.Y] = meshgrid(0:max(ceil(k.range(1)*k.wradius/grid{s}.dx),ceil(k.range(2)*k.wradius/grid{s}.dy)));% grid{s} of searching windows
        [k.ss.el.X_T, k.ss.el.Y_T]=rotredtrans(k.ss.el.X*grid{s}.dx, k.ss.el.Y*grid{s}.dy, k.rotation, k.range); % transforms the grid{s}
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
    Y{scale_i}.sim.xy=grid{s}.xy(isnan(Y{scale_i}.m{1}));
    Y{scale_i}.sim.n=numel(Y{scale_i}.sim.xy);
    if parm.cstk
        Y{scale_i}.sim.xy_r=Y{scale_i}.sim.xy(randperm(Y{scale_i}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
        [Y{scale_i}.sim.x_r,Y{scale_i}.sim.y_r] = ind2sub([grid{s}.ny, grid{s}.nx],Y{scale_i}.sim.xy_r); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
    else
        for i_realisation=1:parm.n_realisation
            Y{scale_i}.sim.xy_r{i_realisation}=Y{scale_i}.sim.xy(randperm(Y{scale_i}.sim.n)); % randomly permutate the ordered vector of index of Y.xy
            [Y{scale_i}.sim.x_r{i_realisation}, Y{scale_i}.sim.y_r{i_realisation}] = ind2sub([grid{s}.ny, grid{s}.nx],Y{scale_i}.sim.xy_r{i_realisation}); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
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
            Y{scale_i}.pt.y = Y{scale_i}.sim.x_r(i_pt); % Find the current point position on the grid{s}. It will be used on each realisation.
            Y{scale_i}.pt.x = Y{scale_i}.sim.y_r(i_pt);
            k0 = kringing_coef(Y{scale_i},X,k,parm,1); % Kriging.
            % variance of the kriging
            Y{scale_i}.pt.s = k0.s;
        end
        
        % * *LIKELIHOOD*
        % We first find the secondary value (and error) (Z.d, Z.std) to
        % create a pdf. This pdf is then multiply inside the kernel to get
        % the density. The likelihood is only...
        Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x))  , Z.std(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x)));
        Y{scale_i}.pt.dens = bsxfun(@times, kernel.dens, Z.pt.dist./sum(Z.pt.dist));
        Y{scale_i}.pt.likelihood = sum(Y{scale_i}.pt.dens, 2)/sum(Y{scale_i}.pt.dens(:));

        for i_realisation=1:parm.n_realisation

            %%
            % * *KRIGING*        
            if ~parm.cstk
                Y{scale_i}.pt.y = Y{scale_i}.sim.x_r{i_realisation}(i_pt); % Find the current point position on the grid{s}. changing for each realisation
                Y{scale_i}.pt.x = Y{scale_i}.sim.y_r{i_realisation}(i_pt);
                k0 = kringing_coef(Y{scale_i},X,k,parm,i_realisation);
                % variance of the kriging
                Y{scale_i}.pt.s = k0.s;
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


            %%
            % * * N-SCORE PRIOR DISTRIBUTION*
            % Back transform the normal distribution (from krigeage) in original space in the kernel.y grid{s}. this become the prior distribution
            Y{scale_i}.pt.prior = Nscore.dist(Y{scale_i}.pt.m, sqrt(Y{scale_i}.pt.s));

            
            %%
            % * *POSTERIORI*
            % Multiply the (nomalized) prior and likelihood to get the (normalised) posteriori
            Y{scale_i}.pt.post_pdf = Y{scale_i}.pt.prior;%.* Y{scale_i}.pt.likelihood;
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
            if parm.plot.krig && i_realisation==1 && ( mod(i_plot+999,1000)==0  || Nscore.forward(Y{scale_i}.pt.sampled)<-3  || Nscore.forward(Y{scale_i}.pt.sampled)>3 ) % 
                figure(1); clf

                subplot(3,2,[1 4])
                imagesc(Y{scale_i}.x,Y{scale_i}.y,Y{scale_i}.m_ns{1}); axis tight; lim=axis;hold on

                tt=-pi:0.01:pi;
                x=Y{scale_i}.x(Y{scale_i}.pt.x)+k.range(1)*cos(tt);
                x2=Y{scale_i}.x(Y{scale_i}.pt.x)+k.wradius*k.range(1)*cos(tt);
                y=Y{scale_i}.y(Y{scale_i}.pt.y)+k.range(2)*sin(tt);
                y2=Y{scale_i}.y(Y{scale_i}.pt.y)+k.wradius*k.range(2)*sin(tt);
                plot(x,y,'--r'); plot(x2,y2,'.-r'); axis(lim)

                lambda_c= 36+60.*(abs(k0.lambda)-min(abs(k0.lambda)))./range(abs(k0.lambda));

                plot(X.x,X.y,'x')
                
                if parm.cstk
                    sel_g=[X.x(k0.sb_mask) X.y(k0.sb_mask); Y{scale_i}.X(k0.ss_mask) Y{scale_i}.Y(k0.ss_mask)];
                    XY_ns = [X.d_ns(k0.sb_mask) ; Y{scale_i}.m_ns{i_realisation}(k0.ss_mask)];
                    scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns,'o','filled','MarkerEdgeColor','w');
                    scatter(Y{scale_i}.x(Y{scale_i}.pt.x),Y{scale_i}.y(Y{scale_i}.pt.y),100,k0.lambda'* XY_ns,'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
                else
                    XY_ns = [X.d_ns; Y{scale_i}.m_ns{1}(~isnan(Y{scale_i}.m_ns{1}))];
                    sel_g_ini=[X.x X.y; Y{scale_i}.X(~isnan(Y{scale_i}.m{i_realisation})) Y{scale_i}.Y(~isnan(Y{scale_i}.m{i_realisation}))];
                    sel_g = sel_g_ini(k0.mask,:);
                    scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns(k0.mask),'o','filled','MarkerEdgeColor','k');
                    scatter(Y{scale_i}.x(Y{scale_i}.pt.x),Y{scale_i}.y(Y{scale_i}.pt.y),100,k0.lambda'* XY_ns(k0.mask),'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
                end
                xlabel('x[m]');ylabel('y[m]');colorbar;

                subplot(3,2,5); hold on;
                plot( kernel.y,Y{scale_i}.pt.prior)
                plot( kernel.y,Y{scale_i}.pt.likelihood)
                plot( kernel.y,Y{scale_i}.pt.post_pdf)
                plot(  kernel.y, max(Y{scale_i}.pt.post_pdf)*Y{scale_i}.pt.post_cdf)
                plot([Y{scale_i}.pt.sampled Y{scale_i}.pt.sampled],[0 max(Y{scale_i}.pt.post_pdf)],'r')
                legend('Prior','Likelihood','Posteriori')


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
            Y{scale_i}.m_ns{i_realisation}(Y{scale_i}.pt.y,Y{scale_i}.pt.x) = Nscore.forward(Y{scale_i}.pt.sampled); % add only the last simulated point to normal score

        end
    end

   
    % display info
    disp(['Simulation ' num2str(scale_i) '/' num2str(numel(parm.scale)) ' finished in : ' num2str(toc(t.tic.scale))])
    t.scale{scale_i} = toc(t.tic.scale);
end

% save intial value
X = X_ini;
if parm.neigh; k.sb.mask = k.sb.mask_ini; end
clear X_ini k.sb.mask_ini

t.global = toc(t.tic.global);
clear t.tic


%% 
% * *SAVE IT*
if parm.saveit
    save(['result/', parm.name ,'_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mat'], 'parm', 'Y', 'grid', 't', 'X', 'Z', 'X_true', 'k', 'kernel', 'Nscore')
end

end