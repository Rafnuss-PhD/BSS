function [krig_m, krig_s]=kringing(Y,X,k)
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

% Super Grid Block from Hard Data:
k.sb.mask_0=k.sb.mask(min([round((Y.y(Y.pt.y)-k.sb.y(1))/k.sb.dy +1)'; k.sb.ny]),  min([round((Y.x(Y.pt.x) -k.sb.x(1))/k.sb.dx +1)'; k.sb.nx])   , :);
sel_sb=[X.x(k.sb.mask_0) X.y(k.sb.mask_0) X.d_ns(k.sb.mask_0)];


% Spiral search per quandrant
[ny,nx]=size(Y.m_ns);
nn_max=length(k.ss.el.dist_s);
n=[0 0 0 0];
sel_ss=cell(4,1);
for q=1:4
    sel_ss{q}=nan(k.nb_neigh(2,q),3);
    nn=2; % 1 is the point itself... therefore unknown
    while n(q)<k.nb_neigh(2,q) && nn<=nn_max && k.ss.el.dist_s(nn)<=k.wradius % while not exceed number of point wanted and still inside the ellipse
        it = Y.pt.x + k.qs(q,1)*k.ss.el.X_s(nn);
        jt = Y.pt.y + k.qs(q,2)*k.ss.el.Y_s(nn);
        if it>0 && jt>0 && it<=nx && jt <=ny % check to not be outside the grid
            if ~isnan(Y.m_ns(jt,it)) % check if it,jt exist
                n(q)=n(q)+1;
                sel_ss{q}(n(q),:) = [Y.X(jt,it) Y.Y(jt,it) Y.m_ns(jt,it)];
            end
        end
        nn=nn+1;
    end
    sel_ss{q}=sel_ss{q}(1:n(q),:); % only the max number of point found.
end


% Combine SuperBlock Point and Spiral Search point.
sel_g=[sel_sb;unique([sel_ss{1};sel_ss{2};sel_ss{3};sel_ss{4}],'rows')];


n=size(sel_g,1);
if any(n<10)
    error(['Not enough point for kriging: ' num2str(n)])
end

% figure; hold on;
% imagesc(Y.x,Y.y,Y.m_ns);
% plot(X.x,X.y,'x')
% plot(Y.x(Y.pt.x),Y.y(Y.pt.y),'o','linewidth',10)
% plot(sel_g(:,1),sel_g(:,2),'o')
% plot(sel_sb(:,1),sel_sb(:,2),'x','linewidth',3)
% legend('Y.m_ns','pt estimated','hard data selected','all data selected')
axis equal



%%
% * *KRIGING*: Find his kringing value in noraml space:
dist_eucl = @(x,y)  sqrt((x/k.range(1)).^2 + (y/k.range(2)).^2);
a0_dist = dist_eucl(sel_g(:,1)-Y.x(Y.pt.x), sel_g(:,2)-Y.y(Y.pt.y));
ab_dist = dist_eucl(bsxfun(@minus,sel_g(:,1),sel_g(:,1)'),  bsxfun(@minus,sel_g(:,2), sel_g(:,2)') );

a0_C=0*a0_dist;
a0_C(a0_dist<k.wradius)=k.C_fx(a0_dist(a0_dist<k.wradius));
ab_C=0*ab_dist;
ab_C(ab_dist<k.wradius)=k.C_fx(ab_dist(ab_dist<k.wradius));

% % Ordinary"Kringing
% lambda = [ab_C ones(sum(n),1); ones(1,sum(n)) 0] \ [a0_C; 1]; 
% krig_m = lambda(1:end-1)'*sel_g(:,3);
% krig_m = lambda(1:end-1)'*sel_g(:,3);
% krig_s = k.var - lambda(1:end-1)'*a0_C - lambda(end);

% Simple Kriging
lambda = ab_C \ a0_C; % Ordinary
krig_m = lambda'*sel_g(:,3);
krig_s = k.var + k.C_0 - lambda'*a0_C;

assert(~isnan(krig_m),'the kriging result is NaN')
assert(krig_s>0,'the kriging std result is less than zero')
end