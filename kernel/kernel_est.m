function kernel = kernel_est(X,Z)
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
%         - n   : Primary variable (1D-X.n elts)


Zdx=nan(size(X.d));
Zstdx=Zdx;
for i=1:X.n
    Zdx(i)=Z.d(X.y(i)==Z.y, X.x(i)==Z.x);
    Zstdx(i)=Z.std(X.y(i)==Z.y, X.x(i)==Z.x);
end



dX=range(X.d)/4; % ???
dZ=3*max(Zstdx(:));% range(Z)/2; % ???

n = 2^8;

% Compute the range for the grid used in the density
min_ZX=[min(Zdx)-dZ, min(X.d)-dX]; 
max_ZX=[max(Zdx)+dZ, max(X.d)+dX];

% Bandwidth proposed by Foster and Bowman
l_z=std(Zdx)* length(Zdx)^(-1/6);
l_x=std(X.d)* length(X.d)^(-1/6);

[~,kernel.dens,x,y]=kde2d_pr([Zdx X.d],n,min_ZX,max_ZX,l_z,l_x);

kernel.x = x(1,:)'; kernel.y = y(:,1); % convert the meshgrid to simple vector
kernel.dy = (kernel.y(2)-kernel.y(1));
kernel.n = n;
kernel.dens(kernel.dens < (8*10^-4)) = 0; % remove very small point WHY??????
end