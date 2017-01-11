%% NSCORE_PERSO
% This function is computing the Normal zscore transform of the input
% vector the function return two fonction handle : one for the normal transform
% and the other one for the back-transform. The function use the inputed
% vector to create the normal transform. Then using interpolation, the
% function handle are created.
%
% INPUT:
% * X            : Input vector
% * support_dist : Vector of value on which the distribution is built
%
% OUTPUT:
% * NscoreT      :
% * NscoreTinv   :
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 29.01.2015

function Nscore = nscore(X,support_dist,method,extrapolationMethod,plotit)

% Compute the empirical cdf
[ft,xt] = ecdf(X);

% ecdf compute the step function like because it is empirical, and there is
% one step per value of X (the height of the step is constent). The real
% pdf of this variable is something (no exacly, but closer) to a linear
% extrapolation between the center of the step. In order to get these
% point, we average the value of x and f:
x = xt(2:end);
f = (ft(1:end-1)+ft(2:end))/2;
% figure; ecdf(X); hold on; plot(x,f,'o-r')


% Now, our problem is to find a good extrapolation at each tail. This is
% important because of that : the prior pdf is estimate in the normal space and back
% transform to initial space to do the sampling. Then the value is
% transform to to normal space for the next prior estimate (kriging).
%
% After trying this:
% x = [support_dist(1)-1 ; x(1)-f(1)*(x(2)-x(1))/((f(2)-f(1)));  x;  x(end)+(1-f(end))*((x(end)-x(end-1))/(f(end)-f(end-1))) ; support_dist(end)+1];
% f = [0+eps; 0+2*eps; f; 1-2*eps; 1-eps];
% x = [support_dist(1); x; support_dist(end)];
%f = [eps; f; 1-eps];
%
% The solution was to fit a expenential, between the first point of the
% kernel and the first point of the ecdf, then use the point of the kernel
% as ...

warning('Implement Hermite Polynomial: PCHIP (see course of Fontainebleau)')
pause;


dx=min([10; numel(X)]); % number of point to take in the interpolation:
if support_dist(1)>0
    fx= 'power1';
elseif support_dist(1)==0
    support_dist(1)=eps;
    fx= 'power1';
else
    fx='exp1';
end

% Lower Tail
w=ones(dx+1,1); w(2)=1000000; % set a very high weight on the last poin to force to be exaclty there (needed for interpolation later)
options = fitoptions('Weights',w); % trick to force the point at x(1)
f_low=fit([support_dist(1); x(1:dx)],[eps; f(1:dx)], fx,options);
% figure;plot(f_low,[support_dist(1); x(1:dx)],[eps; f(1:dx)])

% Upper Tail
w=ones(dx+1,1); w(end-1)=1000000;
options = fitoptions('Weights',w);
f_high=fit([x(end-dx+1:end); support_dist(end)],1-[f(end-dx+1:end); 1-eps],fx,options);
% figure; plot(f_high,[x(end-dx-1:end); support_dist(end)],1-[f(end-dx-1:end); 1-eps])

% Add point upper and lower tail for the interpolation. The point added are
% at support_dist.
x2 = [support_dist(support_dist<x(1))   ; x ; support_dist(support_dist>x(end))];

f2_b=feval(f_low,support_dist(support_dist<x(1)));
f2_e=1-feval(f_high,support_dist(support_dist>x(end)));


f2 = [ f2_b  ; f ; f2_e ];

% check for monoticity in order to interpolate later
i=0;
while ~all(diff(f2)>0)
    id=find(diff(f2)<0);
    f2(id+1) = f2(id+1) + eps;
    i=i+1;
    if i>10
        error('too many correction...')
    end
    f2(id+[-2:2])
end

% f2 and x2 are the vector on which the interpolation between data and
% their probabilite is made.
% figure; ecdf(X); hold on; plot(x,f,'o'); plot(x2,f2,'-d')

Nscore.T_F = griddedInterpolant(x2,f2,method,extrapolationMethod);
Nscore.Tinv_F = griddedInterpolant(f2,x2,method,extrapolationMethod);

% Function of Normal Score Transform. 
% We need to associate to each probability (from the ecdf) a value in the
% normal space which correspond.
Nscore.inverse = @(y) Nscore.Tinv_F(normcdf(y)); % back-transform a value in normal space by taking the normcdf.
% Nscore.forward = @(y) norminv( min([ 1-eps*ones(numel(y),1) max([ eps*ones(numel(y),1) Nscore.T_F(y)],[],2) ],[],2));
Nscore.forward = @(y) norminv( Nscore.T_F(y) );




if plotit
    figure;
    subplot(1,2,1);hist(X); legend(['\mu=' num2str(mean(X)) ' | \sigma=' num2str(std(X))])
    subplot(1,2,2);hist(Nscore.forward(X)); legend(['\mu=' num2str(mean(Nscore.forward(X))) ' | \sigma=' num2str(std(Nscore.forward(X)))])
    
    figure; hold on;
    plot(support_dist,Nscore.forward(support_dist))
    plot(X,Nscore.forward(X),'o')
    xlabel('Initial Space');ylabel('Normal Space');
    legend('support_dist','Hard data X')
    
    figure; hold on;
    [f,x]=hist(X,20); plot(x,f/trapz(x,f));
    [f,x]=hist(Nscore.inverse(randn(400,1)),20); plot(x,f/trapz(x,f));
    xlabel('x'); ylabel('pdf(x)'); legend('Hard data X', 'Nscore.inverse(randn(1000,1))')
    
    figure; hold on;
    ecdf(X)
    ecdf(Nscore.inverse(randn(400,1)))
    ecdf(Nscore.inverse(randn(10000,1)))
    legend('Hard data X', 'Nscore.inverse(randn(400,1))','Nscore.inverse(randn(10000,1))')
end


% Kriging generate a mean and standard deviation in the std normal space.
% We want to transform this std normal distribution in the original space.
% The final distribution is on the grid of support_dist. Therefore we compute
% the nscore transform of the support_dist cell and compute the probability
% corresponding of the normal distribution generated by kriging (mu, sigma)
kernel_y_ns = Nscore.forward(support_dist);

% Method with normpdf
% kernel_b = (kernel_y_ns(3:end)-kernel_y_ns(1:end-2) )/2;
% kernel_b = [kernel_b(1);kernel_b;kernel_b(end)];
% Nscore.dist = @(mu,sigma) normpdf(kernel_y_ns,mu,sigma).*kernel_b;

% Method with normcdf
kernel_y_ns_mid = ( kernel_y_ns(1:end-1)+kernel_y_ns(2:end) ) /2;
Nscore.dist = @(mu,sigma) [ normcdf(kernel_y_ns_mid,mu,sigma) ; 1] - [0 ; normcdf(kernel_y_ns_mid(),mu,sigma)];




if plotit

    figure;
    x=6.2;
    h1=subplot(1,2,1);hold on
    ecdf(X);xlabel('x-original');ylabel('CDF(x)'); xlim([support_dist(1) support_dist(end)])
    plot(x,Nscore.T_F(x),'.r', 'MarkerSize',40)
    line([x x],[0 Nscore.T_F(x)],'Color','k')
    line([x support_dist(end)],[Nscore.T_F(x) Nscore.T_F(x)],'Color','k')
    ax = gca; ax.XTick = sort([ax.XTick ,x]); h1.Position(3)=h1.Position(3)+0.03;
    
    h2=subplot(1,2,2);hold on
    plot(-5:0.1:5,normcdf(-5:0.1:5))
    plot(Nscore.forward(x), Nscore.T_F(x),'.r', 'MarkerSize',40)
    line([-5 Nscore.forward(x)],[Nscore.T_F(x) Nscore.T_F(x)],'Color','k')
    line([Nscore.forward(x) Nscore.forward(x)],[0 Nscore.T_F(x)],'Color','k')
    xlabel('x-Normal Score Transform');ylabel('Standard Normal CDF(x)'); xlim([-5 5])
    ax = gca; ax.XTick = sort([ax.XTick ,Nscore.forward(x)]);ax.YAxisLocation='right';
    h2.Position(1)=h2.Position(1)-0.03;

    
    %%
    mu= 0.4210;
    sigma=1.0044;
    idx=1:20:kernel.n;
    
    figure;
    
    subplot(2,2,1);hold on
    ecdf(X);xlabel('x');ylabel('cdf(x)'); xlim([support_dist(1) support_dist(end)])
    plot(support_dist(idx),Nscore.T_F(support_dist(idx)),'.r', 'MarkerSize',20)
    for i=idx
        line([support_dist(i) support_dist(i)],[0 Nscore.T_F(support_dist(i))],'Color',[0.4,0.4,0.4])
        line([support_dist(i) support_dist(end)],[Nscore.T_F(support_dist(i)) Nscore.T_F(support_dist(i))],'Color',[0.4,0.4,0.4])
    end
    
    subplot(2,2,2);hold on
    plot(-5:0.1:5,normcdf(-5:0.1:5))
    xlabel('x');ylabel('Std Normal cdf(x)'); xlim([-5 5])
    for i=idx
        line([-5 norminv(Nscore.T_F(support_dist(i)))],[Nscore.T_F(support_dist(i)) Nscore.T_F(support_dist(i))],'Color',[0.4,0.4,0.4])
        line([norminv(Nscore.T_F(support_dist(i))) norminv(Nscore.T_F(support_dist(i)))],[0 Nscore.T_F(support_dist(i))],'Color',[0.4,0.4,0.4])
    end
    plot(norminv(Nscore.T_F(support_dist(idx))),Nscore.T_F(support_dist(idx)),'.r', 'MarkerSize',20)
    subplot(2,2,3);hold on
    hist(X);xlabel('x');ylabel('pdf(x)'); xlim([support_dist(1) support_dist(end)])
    ax = gca;ax.YAxisLocation='right';
    
    subplot(2,2,4); hold on;
    plot(-5:0.1:5,normpdf(-5:0.1:5,mu,sigma));
    plot(Nscore.forward(support_dist(idx)),normpdf(Nscore.forward(support_dist(idx)),mu,sigma),'.r', 'MarkerSize',20)
    for i=idx
        line([Nscore.forward(support_dist(i)) Nscore.forward(support_dist(i))],[0 normpdf(Nscore.forward(support_dist(i)),mu,sigma)],'Color','k')
    end
    
    for i=1:length(idx)-1
        u=linspace(Nscore.forward(support_dist(idx(i))),Nscore.forward(support_dist(idx(i+1))),20);
        area([u  Nscore.forward(support_dist(idx(i+1))) Nscore.forward(support_dist(idx(i))) ],...
            [normpdf(u,mu,sigma)  0 0])
    end
    set(gca,'Ydir','reverse'); xlim([-5 5])
    
    
    
    %% 

end


return
%% NOTE:
% Here are 6 Method to compute the transform and back transform
% A:  input vector
% B: Normal z-score of A
% b: point of a std normal distribution
% a: the back transform of b

hold on;

% Method 1: Inital script from Paolo
B=nscoretool(A);
nt=length(A);
zmin=min(b);
zmax=max(b);
ltail=2;
ltpar=2;
utail=1;
utpar=2;
a=backtrtool_pr(b,nt,A,B,zmin,zmax,ltail,ltpar,utail,utpar);

plot(A,B,'o')



% Method 2: mGstat toolbox
% w1,dmin : Extrapolation options for lower tail. w1=1 -> linear interpolation, w1>1 -> Gradual power interpolation
% w2,dmax : Extrapolation options for lower tail. w1=1 -> linear interpolation, w1<1 -> Gradual power interpolation
% DoPlot : plot
% d_nscore : normal score transform of input data
% o_nscore : normal socre object containing information needed to perform normal score backtransform.
[B,o_nscore]=nscore(A,w1,w2,dmin,dmax);
a=inscore(b,o_nscore);



% Method 3. ECDF
CDF_inv=norminv(ecdf(A));
B = CDF_inv(tiedrank(A));
a=quantile(A,normcdf(b));



% Method 4. TieRank
B = norminv( tiedrank(A)/(numel(A)+1));
a=quantile(A,normcdf(b));



% Method 5. http://ch.mathworks.com/help/stats/examples/nonparametric-estimates-of-cumulative-distribution-functions-and-their-inverses.html#zmw57dd0e1074
[Bi,xi] = ecdf(A);
n=numel(A);
xj = xi(2:end);
Bj = (Bi(1:end-1)+Bi(2:end))/2;
xj = [xj(1)-Bj(1)*(xj(2)-xj(1))/((Bj(2)-Bj(1)));  xj;  xj(n)+(1-Bj(n))*((xj(n)-xj(n-1))/(Bj(n)-Bj(n-1)))];
Bj = [0; Bj; 1];

F    = @(y) norminv(interp1(xj,Bj,y,'linear','extrap'));
Finv = @(u) normcdf(interp1(Bj,xj,u,'linear','extrap'));

B=norminv(F(A));
a=Finv(normcdf(b));




% Method 6. http://ch.mathworks.com/help/stats/examples/nonparametric-estimates-of-cumulative-distribution-functions-and-their-inverses.html#zmw57dd0e1147
Fy = ksdensity(A, A, 'function','cdf', 'width',.35);
Finv = @(u) interp1(Fy,y,u,'linear','extrap');

B=norminv(Fy);
a=Finv(normcdf(b));

end