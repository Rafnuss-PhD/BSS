function [k] = covardm_perso(x,x0,covar)

% calculer les covariances
k=zeros(size(x,1),size(x0,1));

for i=1:numel(covar)
    if size(x0,1)==1
        h=sqrt(sum(bsxfun(@minus,x*covar(i).cx,x0*covar(i).cx).^2,2));
        % compute histogram
%         if i==1
%             h_hist = h(:);
%         end
    elseif all(size(x)==size(x0)) && all(x(:)==x0(:))
        h=squareform(pdist(x*covar(i).cx));
    else
        warning('should be here...')
        h=pdist2(x*covar(i).cx, x0*covar(i).cx);
    end

    k=k+kron(covar(i).g(h),covar(i).c0);
end
