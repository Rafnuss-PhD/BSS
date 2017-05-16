addpath(genpath('./.'));
x0=[10 100];
[x fval history] = myproblem(x0);
save('result-BSS/Cal00.mat','history','x','fval');


% parm.par=1;
% parm.parn_n=4;
% parm.n_realisation  = parm.par_n*2;
% parm.notify = 0;
% parm.aggr.method='rad';
% 
%   
% 
% fun = @(x) OF_fx(Klog,Sigma,grid_gen,parm, x(1), x(2) );
% 
% 
% fun([10,100])
% 
% %Single min
% options = optimset('MaxFunEvals',2,'OutputFcn',@outfun,'Display','iter');
% 
% [x,fval,exitflag,output] = fminsearch(fun,[10 100],options);
% 
% save('result-BSS/Cal01.mat','parm','x','fval','exitflag','output');


% Many min
% options = optimset('Display','iter');
% [x,fval,exitflag,output] = fminsearch(fun,[.1 10],options)


% % Pareto Front
% nf = 2; % number of objective functions
% N = 5; % number of points for plotting
% x = zeros(N+1,1);
% f = zeros(N+1,nf);
% x0 = 0.5;
% goal=[0.01 .026];
% options = optimoptions('fgoalattain','Display','iter-detailed','MaxFunEvals',3);
% for r = 0:N
%     t = r/N; % 0 through 1
%     weight = [t,1-t];
%     if r~=0
%         x0=x(r,:);
%     end
%     [x(r+1,:),f(r+1,:),attainfactor,exitflag,output,lambda] = fgoalattain(fun,x0,goal,weight,[],[],[],[],[],[],[],options)
% end
