function kernel = kernel_est(X, Z, range, plotit)
% Nonparametric computes the non-parametric density for two vectors of data
% point X and Z.
% INPUT :
%         - X   : First variable
%         - Z   : Second variable
% 
% OUTPUT :
%         - brandwidth   : Primary variable (1D-X.n elts)
%         - density   : Primary variable (1D-X.n elts)
%         - x   : Primary variable (1D-X.n elts)
%         - y   : Primary variable (1D-X.n elts)
%         - n   : Number of ?? (1D-X.n elts)


Zdx=nan(size(X.d));
Zstdx=Zdx;
for i=1:X.n
    Zdx(i)=Z.d(X.y(i)==Z.y, X.x(i)==Z.x);
    Zstdx(i)=Z.std(X.y(i)==Z.y, X.x(i)==Z.x);
end


n = 2^8;

% Compute the range for the grid used in the density
% dX=range(X.d(:)); % ???
% dZ=range(Zstdx(:));% range(Z)/2; % ???
% min_ZX=[min(Zdx)-dZ, min(X.d)-dX/3]; 
% max_ZX=[max(Zdx)+dZ, max(X.d)+dX/3];
assert(min(Zdx)>range.min(1))
assert(min(X.d)>range.min(2))
assert(max(Zdx)<range.max(1))
assert(max(X.d)<range.max(2))
min_ZX = range.min;
max_ZX = range.max;

% Bandwidth proposed by Foster and Bowman
l_z=std(Zdx)* length(Zdx)^(-1/6);
l_x=std(X.d)* length(X.d)^(-1/6);

[~,kernel.dens,x,y]=kde2d_pr([Zdx X.d],n,min_ZX,max_ZX,l_z,l_x);

kernel.x = x(1,:)'; kernel.y = y(:,1); % convert the meshgrid to simple vector
kernel.dy = (kernel.y(2)-kernel.y(1));
kernel.n = n;

% Normalizing... plus cheating...
kernel.dens(kernel.dens<0) = 0;
kernel.dens = kernel.dens./sum(kernel.dens(:));
% kernel.dens(kernel.dens < (8*10^-4)) = 0; % remove very small point WHY??????

if plotit
    figure; hold on;
    imagesc(kernel.x, kernel.y, kernel.dens)
    scatter(Zdx, X.d,[],X.y, 'x');
    xlabel('Secondary variable X'); ylabel('Primary variable Y');
    axis tight; colorbar('Northoutside')
    plot([0 20],[0 20],'--k')
    % set(gca,'xticklabel',{''});  xlabel(''); 
    % set(gca,'yticklabel',{''});  ylabel('');
    keyboard
end


assert(sum(kernel.dens(:))>0,'Error in the kernel esimator')

end