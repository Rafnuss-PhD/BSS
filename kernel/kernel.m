function kern = kernel(Prim, Sec, krange, plotit)
% kernel computes the non-parametric density
% INPUT :
%         - Prim    : Primary variable
%         - Sec     : Secondary variable
%         - range   : range of the joint pdf
%         - plotit  : Bolean value to plot or not
% 
% OUTPUT :
%         - brandwidth   : Primary variable (1D-X.n elts)
%         - density   : Primary variable (1D-X.n elts)
%         - x   : Primary variable (1D-X.n elts)
%         - y   : Primary variable (1D-X.n elts)
%         - n   : Number of ?? (1D-X.n elts)
if (nargin == 2)
    plotit=1;
    krange={};
    krange.min = [min(Sec.d(:))-.2*range(Sec.d(:)) min(Prim.d(:))-.2*range(Prim.d(:))];
    krange.max = [max(Sec.d(:))+.2*range(Sec.d(:)) max(Prim.d(:))+.2*range(Prim.d(:))];
end


% Find the colocated value of Prim in Sec.
if all(Prim.y(:)==Sec.y) && all(Prim.x(:)==Sec.x)
    data=[Sec.d Prim.d];
else
    data = nan(Prim.n,2);
    for i=1:Prim.n
        [~,idy] = min((Prim.y(i)-Sec.y).^2);
        [~,idx] = min((Prim.x(i)-Sec.x).^2);
        data(i,:)=[Sec.d(idy, idx) Prim.d(i)];
    end
end

% Number of point in the joint pdf
n = 2^8;


% Define the range of the grid used in the density. Range have been
% computed as input data
assert( all(max(data) < krange.max) )
assert( all(min(data) > krange.min) )


% Bandwidth proposed by Foster and Bowman
l_sec=std(data(:,1))* Prim.n^(-1/6);
l_prim=std(data(:,2))* Prim.n^(-1/6);

% Compute the density
[~,kern.dens,x,y]=kde2d_pr(data,n,krange.min,krange.max,l_sec,l_prim);

% recover and re-name the output
kern.axis_sec = x(1,:)'; 
kern.axis_prim = y(:,1); % convert the meshgrid to simple vector
kern.daxis_prim = (kern.axis_prim(2)-kern.axis_prim(1));
kern.n = n;

% Normalizing...
kern.dens(kern.dens<0) = 0;
kern.dens = kern.dens./sum(kern.dens(:));

% Compute the marginal distribution
kern.prior = hist(Prim.d,kern.axis_prim)';
kern.prior = kern.prior./sum(kern.prior(:));

% Plot
if plotit
    figure; hold on;
    imagesc(kern.axis_sec, kern.axis_prim, kern.dens)
    scatter(data(:,1), data(:,2),'xk');
    ylabel('Primary variable'); xlabel('Secondary variable');
    axis tight; colorbar('Northoutside')
    % plot([0 20],[0 20],'--k')
    % set(gca,'xticklabel',{''});  xlabel(''); 
    % set(gca,'yticklabel',{''});  ylabel('');
    keyboard
end

% Check for error
assert(sum(kern.dens(:))>0,'Error in the kernel esimator')

end
