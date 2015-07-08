%% BSGS_ONE_SIMULATION is computing one simulation of BSGS.
% It's a sequential stochatic simulation of a primary variable from
% combination of sparce primary (X) data and more dence secondary data (Z).
% The grid is define by the size of the secondary variable (Z).
%
% INPUT:
%
% * X.d     : Primary variable (1D-X.n elts)
% * X.x     : x-coordinate of the primary variable (1D-X.n elts)
% * X.y     : y-coordinatesize(Z.d) of the primary variable (1D-X.n elts)
% * X.n     : number of elements
% * X.xy    : x and y-coordinate in matrice index from grid matric (1D-X.n elts)
% * Z.d     : Secondary variable (2D- grid.nx x grid.n elts)
% * Z.std   : Secondary variable std error (2D- grid.nx x grid.n elts)
% * Y.x     : x-coordinate of the resulting primary variable in the order of simulation
% * Y.y     : y-coordinate of the resulting primary variable in the order of simulation
% * U       : randomly generated value between 0 and 1 to be use to sample in the posteriori distribution
% * kernel  : kernel density information
% * k       : kriging information
%
% OUTPUT:
%
% * Y_mat   : Simulated Primary variable in matrix (3rd dim for simulation)
% * t       : time of simulation
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

function [Y_m, dy] = BSGS_pt(X, Z, Y, U, kernel, Nscore, k, grid)

%%
% * *RESULT ALLOCATION*
n_sim=size(Y.m,1);
Y.m_ns = repmat({nan(size(Y.m{1}))},n_sim,1); % matrix des resutlats
% Maybe this is not very efficient, we should use normal transform for all
% data and only back-tranform at the end of the simulation those of
% interest
I=~isnan(Y.m{1});
for i_sim=1:n_sim
    Y.m_ns{i_sim}(I) = Nscore.forward(Y.m{i_sim}(I));
end
X.d_ns = Nscore.forward(X.d);


%%
% * *OTHER*
ptmod=1000; % display information and check for nscore tool every ptmod point


%%
% * *SEARCHING-WINDOWS*: create the windows for kringing with the function to compute
% the normalized distence and the order of visit of the cells.
% Spiral Search setting: previously data (on grid location)
[k.ss.el.X, k.ss.el.Y] = meshgrid(0:max(ceil(k.range(1)*k.wradius/grid.dx),ceil(k.range(2)*k.wradius/grid.dy)));% grid of searching windows
[k.ss.el.X_T, k.ss.el.Y_T]=rotredtrans(k.ss.el.X*grid.dx, k.ss.el.Y*grid.dy, k.rotation, k.range); % transforms the grid
k.ss.el.dist = sqrt(k.ss.el.X_T.^2 + k.ss.el.Y_T.^2);
[k.ss.el.dist_s, k.ss.el.dist_idx] = sort(k.ss.el.dist(:)); % sort according distence.
k.ss.el.X_s=k.ss.el.X(k.ss.el.dist_idx); k.ss.el.Y_s=k.ss.el.Y(k.ss.el.dist_idx);


%% Point simulation
% Here each point is simulated
dy=zeros(Y.sim.n,1); % What is that ?

for i_pt=1:Y.sim.n;
    
    %%
    % * *PTMOD*: Every PTMOD points, it display information and check for NSCORE
    %     if mod(i_pt,ptmod)==0
    %         disp(['indice ',num2str(i_pt),' over ' ,num2str(Y.sim.n), ' (', num2str(i_pt/Y.sim.n*100),'%) ','  time  ',num2str(toc) ,' s (expected: ',num2str(toc*Y.sim.n/i_pt/60), ' min - and remaining: ',num2str(toc*Y.sim.n/i_pt/60-toc/60),')'])
    %         t(u)=toc;u=u+1;
    %     end
    
    %%
    % * *1. POINT*: Select the point
    Y.pt.y = Y.sim.x_r(i_pt); % Find the current point position on the grid
    Y.pt.x = Y.sim.y_r(i_pt);
    
    
    %%
    % * *3. KRIGING-1*: Find the kriging coefficient
    plotit=0;
    k0 = kringing_coef(Y,X,k,plotit);
    clear plotit
    
    
    for i_sim=1:n_sim
        
        %%
        % * *3. KRIGING-2*: Apply the coefficient
        Y.pt.m = k0.lambda'* [X.d_ns(k0.sb_mask) ; Y.m_ns{i_sim}(k0.ss_mask)];
        Y.pt.s = k0.s;
        
        %%
        % * *4. N-SCORE PRIOR DISTRIBUTION*: Back transform the normal distribution (from krigeage) in original space in the kernel.y grid. this become the prior distribution
        Y.pt.prior = Nscore.dist(Y.pt.m, sqrt(Y.pt.s));
        
        %     if (.99>=sum(Y.pt.prior)) || (sum(Y.pt.prior)>=1.01)
        %         error(['Distribution of the prior is not exacly 1...: ' num2str(sum(Y.pt.prior)) ])
        %         plot(kernel.y, Y.pt.prior)
        %         Y.pt.prior = (Y.pt.prior./sum(Y.pt.prior));
        %     end
        
        
        %%
        % * *6. LIKELIHOOD*:
        
        Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==Y.Y(Y.pt.y,Y.pt.x), Z.x==Y.X(Y.pt.y,Y.pt.x))  , Z.std(Z.y==Y.Y(Y.pt.y,Y.pt.x), Z.x==Y.X(Y.pt.y,Y.pt.x)));
        %     if (sum(Z.pt.dist)-sum(Z.pt.dist(2:end-1)))/sum(Z.pt.dist)>0.01
        %         warning('kernel.x seems to be too small to encompass Z.d at Y.pt distribution')
        %         plot(kernel.x, Z.pt.dist)
        %     end
        Y.pt.dens = bsxfun(@times, kernel.dens, Z.pt.dist./sum(Z.pt.dist));
        Y.pt.likelihood=sum(Y.pt.dens, 2);
        
        %%
        % * *7. POSTERIORI*: Multiply the (nomalized) prior and likelihood
        % to get the (normalised) posteriori
        Y.pt.post_pdf = Y.pt.prior .*(Y.pt.likelihood./sum(Y.pt.likelihood));
        Y.pt.post_pdf = Y.pt.post_pdf./sum(Y.pt.post_pdf);
        
        %%
        % * *8. SAMPLING*: Sample a point in the posteri distribution. To
        % do so the CDF is created and interpolation tool is used to assure
        % randomized sampling. The second part is added to make it mono... increasing for the interpolation to work
        % sort value from cdf of the posterior
        Y.pt.post_cdf = cumsum(Y.pt.post_pdf) +  ((1:1:length(kernel.y))*(1e-10))';
        Y.pt.sampled = interp1(Y.pt.post_cdf, kernel.y, min(Y.pt.post_cdf)+range(Y.pt.post_cdf)*U(i_pt,i_sim));
        Y.m{i_sim}(Y.pt.y,Y.pt.x) = Y.pt.sampled;
        Y.m_ns{i_sim}(Y.pt.y,Y.pt.x) = Nscore.forward(Y.pt.sampled); % add only the last simulated point to normal score

        
    end
end
% save('M')
Y_m=Y.m;
end