function [Y,t_scale]=BSGS_par_in(X, Z, kernel, k, Nscore, grid, parm)

Y=cell(parm.n_scale,1); % allocate variable space

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
                % We first find the secondary value (and error) (Z.d, Z.std) to create a pdf. This pdf is then multiply inside the kernel to get
                % the density. The likelihood is only...
                Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x))  , Z.std(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x)));
                Y{scale_i}.pt.dens = bsxfun(@times, kernel.dens, Z.pt.dist./sum(Z.pt.dist));
                Y{scale_i}.pt.likelihood = sum(Y{scale_i}.pt.dens, 2)/sum(Y{scale_i}.pt.dens(:));
                Z.pt.dist = normpdf(kernel.x, Z.d(Z.y==,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x))  , Z.std(Z.y==Y{scale_i}.Y(Y{scale_i}.pt.y,Y{scale_i}.pt.x), Z.x==Y{scale_i}.X(Y{scale_i}.pt.y,Y{scale_i}.pt.x)));
                
            end
            
            if parm.neigh % if option smart neihbour is selected
                % the weight(lambda) were computed earlier in kriging_coef.m,
                % then the point are coming from the hard data (X.d_ns) and the
                % previously simulated point Y{scale_i}.m_ns.
                Y{scale_i}.pt.m = k0.lambda'* [X.d_ns(k0.sb_mask) ; Y{scale_i}.m_ns{i_realisation}(k0.ss_mask)];
                if isempty(Y{scale_i}.pt.m)
                    warning('Empty estimate. replace by 0')
                    Y{scale_i}.pt.m=0;
                end
            else
                XY_ns = [X.d_ns; Y{scale_i}.m_ns{i_realisation}(~isnan(Y{scale_i}.m_ns{i_realisation}))];
                Y{scale_i}.pt.m = k0.lambda'* XY_ns(k0.mask);
            end
            
            
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
            Y{scale_i}.pt.post_pdf =  Y{scale_i}.pt.likelihood.^parm.w_X(scale_i) .* Y{scale_i}.pt.prior.^parm.p_w(2,scale_i);
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
            if scale_i==9 && parm.plot.krig && i_realisation==1 && (i_plot==1 || mod(i_plot+49,50)==0  || Nscore.forward(Y{scale_i}.pt.sampled)<-3  || Nscore.forward(Y{scale_i}.pt.sampled)>3 ) % 
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
                %[f,x]=hist(X_ini.d(:)); plot(x,f/trapz(x,f),'linewidth',2);
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
    disp(['Simulation ' num2str(scale_i) '/' num2str(parm.n_scale) ' finished in : ' num2str(toc(t.tic.scale))])
    t_scale{scale_i} = toc(t.tic.scale);
end

end
