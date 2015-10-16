%% BSGS_S is computing simulation of BSGS at one specified scale.
% It's a sequential stochatic simulation of a primary variable from
% combination of sparce primary (X) data and more dence secondary data (Z).
% The grid is define by the size of the secondary variable (Z).
%
% INPUT:
%
% * X       : Primary variable (1D-X.n elts)
% * Z       : Secondary variable (2D- grid.nx x grid.n elts)
% * Y       : Simulation (2D- grid.nx x grid.n elts)
% * Y.x     : x-coordinate of the resulting primary variable in the order of simulation
% * Y.y     : y-coordinate of the resulting primary variable in the order of simulation
% * kernel  : kernel density information
% * Nscore  :
% * k       : kriging information
% * grid    : grid information
% * parm    : parameter
%
% OUTPUT:
%
% * Y_m     : Simulated Primary variable
% * t       : time of simulation
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

function [Y_m, Y_m_ns] = BSGS_s(X, Z, Y, kernel, Nscore, k, grid, parm)

%%
% * *RESULT ALLOCATION*
% Create empty matrix for resultats in normal space. Then populate it with
% known value
Y.m_ns = repmat({nan(size(Y.m{1}))},parm.n_realisation,1); 

I=~isnan(Y.m{1});
for i_realisation=1:parm.n_realisation
    Y.m_ns{i_realisation}(I) = Nscore.forward(Y.m{i_realisation}(I));
end

X.d_ns = Nscore.forward(X.d);


%%
% * *SEARCHING-WINDOWS*: create the windows for kringing with the function to compute
% the normalized distence and the order of visit of the cells.
% Spiral Search setting: previously data (on grid location)
if parm.neigh
    [k.ss.el.X, k.ss.el.Y] = meshgrid(0:max(ceil(k.range(1)*k.wradius/grid.dx),ceil(k.range(2)*k.wradius/grid.dy)));% grid of searching windows
    [k.ss.el.X_T, k.ss.el.Y_T]=rotredtrans(k.ss.el.X*grid.dx, k.ss.el.Y*grid.dy, k.rotation, k.range); % transforms the grid
    k.ss.el.dist = sqrt(k.ss.el.X_T.^2 + k.ss.el.Y_T.^2);
    [k.ss.el.dist_s, k.ss.el.dist_idx] = sort(k.ss.el.dist(:)); % sort according distence.
    k.ss.el.X_s=k.ss.el.X(k.ss.el.dist_idx); k.ss.el.Y_s=k.ss.el.Y(k.ss.el.dist_idx);
end


%% Generate the order for visting cells
% Randomly permute the cell not known (to be visited). And generate the
% Y.x and Y.y coordinate in a random order.
Y.sim.xy=grid.xy(isnan(Y.m{1}));
Y.sim.n=numel(Y.sim.xy);
if parm.cstk
    Y.sim.xy_r=Y.sim.xy(randperm(Y.sim.n)); % randomly permutate the ordered vector of index of Y.xy
    [Y.sim.x_r,Y.sim.y_r] = ind2sub([grid.ny, grid.nx],Y.sim.xy_r); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
else
    for i_realisation=1:parm.n_realisation
        Y.sim.xy_r{i_realisation}=Y.sim.xy(randperm(Y.sim.n)); % randomly permutate the ordered vector of index of Y.xy
        [Y.sim.x_r{i_realisation}, Y.sim.y_r{i_realisation}] = ind2sub([grid.ny, grid.nx],Y.sim.xy_r{i_realisation}); % * ind2sub() is taking linearized matrix index (i) and transform matrix index (i,j). We insert the randomized path of this specific simulation (ns) in Y.x and Y.y after the already known data of the primary data
    end
end

%%
% * *RANDOM NORMAL FIELD* 
% This create the random Normal distribution used for sampling the posteriori distribution at each point.
U = normcdf(randn(Y.sim.n,parm.n_realisation)); % random field


%% Point simulation

i_plot=0; % used for computing the number of point simulated. Used for ploting
for i_pt=1:Y.sim.n; % loop over each point
    
    i_plot=i_plot+1;

    if parm.cstk % if option constant weight is activate.
        Y.pt.y = Y.sim.x_r(i_pt); % Find the current point position on the grid. It will be used on each realisation.
        Y.pt.x = Y.sim.y_r(i_pt);
        k0 = kringing_coef(Y,X,k,parm,1); % Kriging.
    end
    
    for i_realisation=1:parm.n_realisation
        
        if ~parm.cstk
            Y.pt.y = Y.sim.x_r{i_realisation}(i_pt); % Find the current point position on the grid. changing for each realisation
            Y.pt.x = Y.sim.y_r{i_realisation}(i_pt);
            k0 = kringing_coef(Y,X,k,parm,i_realisation);
        end

        if parm.neigh % if option smart neihbour is selected
            % the weight(lambda) were computed earlier in kriging_coef.m,
            % then the point are coming from the hard data (X.d_ns) and the
            % previously simulated point Y.m_ns.
            Y.pt.m = k0.lambda'* [X.d_ns(k0.sb_mask) ; Y.m_ns{i_realisation}(k0.ss_mask)];
        else
            XY_ns = [X.d_ns; Y.m_ns{i_realisation}(~isnan(Y.m_ns{i_realisation}))];
            Y.pt.m = k0.lambda'* XY_ns(k0.mask);
        end
        
        % variance of the kriging
        Y.pt.s = k0.s;
        
        
        %%
        % * *LIKELIHOOD*
        % We first find the secondary value (and error) (Z.d, Z.std) to
        % create a pdf. This pdf is then multiply inside the kernel to get
        % the density. The likelihood is only...
        Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==Y.Y(Y.pt.y,Y.pt.x), Z.x==Y.X(Y.pt.y,Y.pt.x))  , Z.std(Z.y==Y.Y(Y.pt.y,Y.pt.x), Z.x==Y.X(Y.pt.y,Y.pt.x)));
        Y.pt.dens = bsxfun(@times, kernel.dens, Z.pt.dist./sum(Z.pt.dist));
        Y.pt.likelihood = sum(Y.pt.dens, 2)/sum(Y.pt.dens(:));
        
        %%
        % * *NEW GRID p*
        % This need to be transform in normal space to be combine with
        % kriging estimate. We difine a new grid p  in normal space for 
        % both the likelihood (interpolation is needed) and kriging.
        % p is limited at the coordinate which reprensent 0.1% of
        % probabllity
        p.lim=0.001;
        p.min = min([ kernel.y_ns(find(Y.pt.likelihood>p.lim,1,'first'))   ;  norminv(p.lim,Y.pt.m,Y.pt.s)]);
        p.max = max([ kernel.y_ns(find(Y.pt.likelihood>p.lim,1,'last'))   ;  norminv(1-p.lim,Y.pt.m,Y.pt.s)]);
        p.d = linspace(p.min,p.max,500);
        
        Y.pt.likelihood = interp1(kernel.y_ns(kernel.y_ns_unique), Y.pt.likelihood(kernel.y_ns_unique), p.d,'pchip');
        Y.pt.likelihood = Y.pt.likelihood /sum(Y.pt.likelihood);
        Y.pt.prior = normpdf(p.d,Y.pt.m,Y.pt.s);
        Y.pt.prior = Y.pt.prior./sum(Y.pt.prior);

        %%
        % * *7. POSTERIORI*: Multiply the (nomalized) prior and likelihood
        % to get the (normalised) posteriori
        Y.pt.post_pdf = Y.pt.prior;%.* Y.pt.likelihood;
        Y.pt.post_pdf = Y.pt.post_pdf./sum(Y.pt.post_pdf);
        
        %%
        % * *8. SAMPLING*: Sample a point in the posteri distribution. To
        % do so the CDF is created and interpolation tool is used to assure
        % randomized sampling. 
        Y.pt.post_cdf = cumsum(Y.pt.post_pdf);

        % Interpolation only works if it is mono.. increasing. But the cdf
        % can reach 1 so we just add eps (very small value).
        if ~all(diff(Y.pt.post_cdf)>0) 
            Y.pt.post_cdf = Y.pt.post_cdf +  linspace(0,numel(Y.pt.post_cdf)*eps*2,numel(Y.pt.post_cdf));
        end
        
        Y.pt.sampled = interp1(Y.pt.post_cdf, p.d, U(i_pt,i_realisation),'pchip');
        
        
        
        %%
        % * *PLOTIT*
        if parm.plot.krig && i_realisation==1 && ( mod(i_plot+999,1000)==0  || Y.pt.sampled<-3  || Y.pt.sampled>3 ) % 
            figure(1); clf
          
            subplot(3,2,[1 4])
            imagesc(Y.x,Y.y,Y.m_ns{1}); axis tight; lim=axis;hold on
            
            t=-pi:0.01:pi;
            x=Y.x(Y.pt.x)+k.range(1)*cos(t);
            x2=Y.x(Y.pt.x)+k.wradius*k.range(1)*cos(t);
            y=Y.y(Y.pt.y)+k.range(2)*sin(t);
            y2=Y.y(Y.pt.y)+k.wradius*k.range(2)*sin(t);
            plot(x,y,'--r'); plot(x2,y2,'.-r'); axis(lim)
            
            lambda_c= 36+60.*(abs(k0.lambda)-min(abs(k0.lambda)))./range(abs(k0.lambda));

            if parm.cstk
                sel_g=[X.x(k0.sb_mask) X.y(k0.sb_mask); Y.X(k0.ss_mask) Y.Y(k0.ss_mask)];
                XY_ns = [X.d_ns(k0.sb_mask) ; Y.m_ns{i_realisation}(k0.ss_mask)];
                scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns,'o','filled','MarkerEdgeColor','w');
                scatter(Y.x(Y.pt.x),Y.y(Y.pt.y),100,k0.lambda'* XY_ns,'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
            else
                XY_ns = [X.d_ns; Y.m_ns{1}(~isnan(Y.m_ns{1}))];
                sel_g_ini=[X.x X.y; Y.X(~isnan(Y.m{i_realisation})) Y.Y(~isnan(Y.m{i_realisation}))];
                sel_g = sel_g_ini(k0.mask,:);
                scatter(sel_g(:,1),sel_g(:,2),lambda_c,XY_ns(k0.mask),'o','filled','MarkerEdgeColor','k');
                scatter(Y.x(Y.pt.x),Y.y(Y.pt.y),100,k0.lambda'* XY_ns(k0.mask),'o','filled','MarkerEdgeColor','r','LineWidth',1.5)
            end
            xlabel('x[m]');ylabel('y[m]');colorbar;
            
            subplot(3,2,5); hold on;
            plot(p.d,Y.pt.prior)
            plot(p.d,Y.pt.likelihood)
            plot(p.d,Y.pt.post_pdf)
            plot( p.d, max(Y.pt.post_pdf)*Y.pt.post_cdf)
            plot([Y.pt.sampled Y.pt.sampled],[0 max(Y.pt.post_pdf)],'r')
            legend('Prior','Likelihood','Posteriori')
            
            
            subplot(3,2,6); hold on;
            [f,x]=hist(X.d_ini(:),20); plot(x,f/trapz(x,f),'linewidth',2);
            [f,x]=hist(Z.d(:),40); plot(x,f/trapz(x,f),'linewidth',2);
            [f,x]=hist(Y.m{i_realisation}(:),50); plot(x,f/trapz(x,f));
            legend('Well sampling', 'ERT','Simulation')
            
            drawnow
            keyboard
        end
        
        %% Back transform 
        % * *
        Y.m{i_realisation}(Y.pt.y,Y.pt.x) = Nscore.inverse(Y.pt.sampled);
        Y.m_ns{i_realisation}(Y.pt.y,Y.pt.x) = Y.pt.sampled; % add only the last simulated point to normal score
        
    end
end

Y_m = Y.m;
Y_m_ns = Y.m_ns;
end