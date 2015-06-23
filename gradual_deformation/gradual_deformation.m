%% Gradual Deformation
% In this case, we use a very simple graduation simulation. We start with
% the first fiel as the best field. Then we deforme the best field with the
% challenger. The result from the deformation become the best field.
%
% INPUT:
%
% * X       : Primary variable
% * Z       : Secondary variable
% * Y_true  : True value of Y
% * U       : randomly generated value between 0 and 1 to be use to sample in the posteriori distribution
% * kernel  : kernel density information
% * kri     : kriging information
%
% OUTPUT:
%
% * Y_mat   : Simulated Primary variable in matrix (3rd dim for simulation)
% * t_sim   : time of simulation
%
% * *Author:* Raphael Nussbaumer (raphael.nussbaumer@unil.ch)
% * *Date:* 02.02.2015

function [Y_mat_best, U_best] = gradual_deformation(X, Z, Y_sim, U,  kernel, k, Y_true)
%%
% * *INFO*: get some info
[Y.n_sim, n_sim] = size(U);

%%
% * *FX HANDLE* : Build some function handle: OF_fx (objectif fonction
% (RMSE)), Grad_fx (gradual deformation), OF_Grad_fx (Objetif function gloabal)
OF_fx = @(Y) sqrt(mean(mean((Y_true-Y).^2))); 
Grad_fx = @(U1,U2,theta) U1*cos(theta) + U2*sin(theta);
OF_Grad_fx = @(U1,U2,theta) OF_fx(BSGS_one_simulation(X, Z, Y_sim, normcdf(Grad_fx(U1,U2,theta)), kernel, k));

%%
% * *OPTIMIZATION ALLOCATION*
U_best = zeros(Y.n_sim,n_sim); % store the best U at each iteration
U_best(:,1) = U(:,1);
OF_best = zeros(n_sim,1); % store the best OF at each iteration
OF_best(1) = OF_fx(BSGS_one_simulation(X, Z, Y_sim, normcdf(U_best(:,1)), kernel, k));

%%
% * *MINIMIZE FUNCTION*: option selection
options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',5,'PlotFcns',@optimplotfval,'TolX',1e-12);

%%
% * *ITERATIVE GRADUAL DEFORMATION*
for ns=2:n_sim
    fun_min = @(theta) OF_Grad_fx(U_best(:,ns), U(:,ns), theta); % Fun to minimized
    [theta_best,OF_best_fminbnd] = fminbnd(fun_min,0,pi/2,options); %matlab function to find the minimal
    
    if OF_best_fminbnd < OF_best(ns-1) % if fminbnd improve the OF score, take this new score
        OF_best(ns) = OF_best_fminbnd;
        U_best(:,ns) = Grad_fx(U_best(:,1), U(:,ns), theta_best);
    else % otherwise keep the previous
        OF_best(ns) = OF_best(ns-1);
        U_best(:,ns) = U_best(:,ns-1);
    end
end

Y_mat_best = BSGS_one_simulation(X, Z, Y_sim, normcdf(U_best(:,ns)), kernel, k);

stop
end