function k0=kringing_coef(Y,X,k,plotit)
%% kriging return the mean and variance estimate at the point Y.pt
%
% INPUT:
%
% * Y       : Primary variable
% * X       : Primary variable
% * k       : kriging information
%
% OUTPUT:
%
% * krig_m  : kriging mean estimate
% * krig_s  : kriging variance estimate
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015


%% *SELECTION OF DATA*
% Use Supergrid Block for hard data and spiral search for previously
% simulated point. combine all the point and check few stuff

% Super Grid Block from Hard Data:
k0.sb_mask=reshape(k.sb.mask(min([round((Y.y(Y.pt.y)-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny]),  min([round((Y.x(Y.pt.x) -k.sb.x(1))/k.sb.dx +1)'; k.sb.nx])   , :),X.n,1);


% Spiral search per quandrant
nn_max=length(k.ss.el.dist_s);
n=[0 0 0 0];
sel_ss_idx=cell(4,1);

for q=1:4
    sel_ss_idx{q}=nan(k.nb_neigh(2,q),2);
    nn=2; % 1 is the point itself... therefore unknown
    while n(q)<k.nb_neigh(2,q) && nn<=nn_max && k.ss.el.dist_s(nn)<=k.wradius % while not exceed number of point wanted and still inside the ellipse
        it = Y.pt.x + k.qs(q,1)*k.ss.el.X_s(nn);
        jt = Y.pt.y + k.qs(q,2)*k.ss.el.Y_s(nn);
        if it>0 && jt>0 && it<=Y.nx && jt <=Y.ny % check to not be outside the grid
            if ~isnan(Y.m_ns{1}(jt,it)) % check if it,jt exist
                n(q)=n(q)+1;
                sel_ss_idx{q}(n(q),:) = [jt it];
            end
        end
        nn=nn+1;
    end
    sel_ss_idx{q}=sel_ss_idx{q}(1:n(q),:); % only the max number of point found.
end

k0_ss_idx = unique([sel_ss_idx{1};sel_ss_idx{2};sel_ss_idx{3};sel_ss_idx{4}],'rows');
k0.ss_mask = sub2ind([Y.ny, Y.nx],k0_ss_idx(:,1),k0_ss_idx(:,2));


% Combine SuperBlock Point and Spiral Search point.
sel_g=[X.x(k0.sb_mask) X.y(k0.sb_mask); Y.X(k0.ss_mask) Y.Y(k0.ss_mask)];

if ~(size(unique(sel_g,'rows'),1)==size(sel_g,1))
    error('None unique point for kriging: ')
end


if size(sel_g,1)<10
    error(['Not enough point for kriging: ' num2str(n)])
end



%%
% * *KRIGING*: Find his kringing value in noraml space:


a0_C=covardm(sel_g,[Y.x(Y.pt.x) Y.y(Y.pt.y)],k.model,k.var);
ab_C=covardm(sel_g,sel_g,k.model,k.var);


k0.lambda = ab_C \ a0_C; % Ordinary
k0.s = sum(k.var) - k0.lambda'*a0_C;

assert(any(~isnan(k0.lambda)),'the kriging coeff is NaN')
assert(k0.s>0,'the kriging std result is less than zero')

%%
% * *PLOTIT*
if plotit
    figure('units','normalized','outerposition',[0 0 1 1]); hold on
    imagesc(Y.x,Y.y,Y.m_ns{1});
    plot(Y.x(Y.pt.x),Y.y(Y.pt.y),'o','linewidth',10)
    lambda_c= min(Y.m_ns{1}(:))   +  range(Y.m_ns{1}(:)).*(k0.lambda-min(k0.lambda))./range(k0.lambda);
    scatter(X.x(k0.sb_mask), X.y(k0.sb_mask),[],lambda_c(1:sum(k0.sb_mask)),'o','filled')
    scatter(Y.X(k0.ss_mask), Y.Y(k0.ss_mask),[],lambda_c(sum(k0.sb_mask)+1:end),'d','filled')
    legend('pt estimated','hard data selected','soft data selected')
    axis tight
    pause
end

end