
fun = @(x,mu,sigma) 1/sqrt(2*pi*sigma^2) * exp(-(x-mu).^2/(2*sigma^2));

mu=[1 1];
s=[1 1];
w=[1 1];

fun2 = @(x) fun(x,mu(1),s(1)).^w(1)  .* fun(x,mu(2),s(2)).^w(2);
fun3 = @(u) integral(fun2,-Inf,u) / integral(fun2,-Inf,Inf);


funA = @(u) integral( @(x) (1/sqrt(2*pi*s(1)^2) * exp(-(x-mu(1)).^2/(2*s(1)^2))).^w(1)  .* (1/sqrt(2*pi*s(2)^2) * exp(-(x-mu(2)).^2/(2*s(2)^2))).^w(2),-Inf,u,'RelTol',0.1,'ArrayValued',true);


opt = optimset('Display','off','TolFun',.01);
n=(513*65)*1;
u=nan(1,n);
tic
for i=1:n
    u(i) = fsolve(funA,rand(1),opt);
end
toc

xx=-3:0.1:3;

figure; hold on;
plot(xx,fun(xx,mu(1),s(1)))
plot(xx,fun(xx,mu(2),s(2)))
plot(xx,fun2(xx))



syms x

for i=1:n
mu=[1 1];
s=[1 1];
w=[1 1];

fun = (1/sqrt(2*pi*s(1)^2) * exp(-(x-mu(1)).^2/(2*s(1)^2))).^w(1)  .* (1/sqrt(2*pi*s(2)^2) * exp(-(x-mu(2)).^2/(2*s(2)^2))).^w(2);

intfun = int(fun,x);

% fplot(intfun);

solve(intfun == 0);
end