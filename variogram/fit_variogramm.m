function k = fit_variogramm(X,Z,parm,X_true)




%% Find the range form Z

[gamma_x, gamma_y] = variogram_gridded_perso(Z.d);
figure; hold on;
plot(Z.x(1:end/2),gamma_x(1:end/2))
plot(Z.y(1:end/2),gamma_y(1:end/2))
id = grid_gen.x<parm.k.range(1)*3;


[gamma_x, gamma_y] = variogram_gridded_perso(phi_true);
figure; hold on;
plot(grid_gen.x(1:end/2),gamma_x(1:end/2))
plot(grid_gen.y(1:end/2),gamma_y(1:end/2))


%% vertical
% LSCurve global fit
myfun = @(x,h) semivariogram1D(h,x(1),x(2),'sph', x(3)); % x= [ c0(1) range c0(2) (nugguet effect)] 
options = optimoptions('lsqcurvefit','Display','off');

% Compute the empirical VERTICAL variogram with sampled point at well
value_y = (X.d(:)-mean(X.d(:)))/std(X.d(:));
coord_y = X.y(:)+1000000*X.x(:); % we add an offset between each line of the 2D data so that each line are corolated separately
dy_avg = mean(diff(sort(unique(X.y(:))))); % compute the average delta y value of the data to use bins of this width
nrbins = ( max(X.y(:))-min(X.y(:)) ) / (dy_avg*2);
Emp = variogram(coord_y,value_y,'nrbins',nrbins,'plotit',false,'maxdist',max(X.y(:)),'subsample',20000);
assert(all(~isnan(Emp.val)),'problem');


range0 = 10;
c0 = [.95 .05];
x0 = [c0(1) range0 c0(2)];
x = lsqcurvefit(myfun,x0,Emp.distance,Emp.val,[0 0 0],[1 5 .05],options);

y_range=x(2);
c=[x(1); x(3)];

% plotit
if plotit
    figure;
    subplot(1,2,1);hold on;
    plot(Emp.distance,Emp.val,'--o')
    plot(Emp.distance,myfun([c(1) y_range c(2)],Emp.distance))
    legend('Empirical Variogram vertical','Fitted Variogram horizontal')
end

%% Compute the HORIZONTAL variogram
value_x = (Z.d_raw(:)-mean(Z.d_raw(:)))/std(Z.d_raw(:));
Z.X_raw = repmat(Z.x_raw',numel(Z.y_raw),1);
Z.Y_raw = repmat(Z.y_raw,1,numel(Z.x_raw));
coord_x = Z.X_raw(:)+1000000*Z.Y_raw(:); % we add an offset between each line of the 2D data so that each line are corolated separately
dx_avg = mean(diff(sort(unique(Z.x_raw(:))))); % compute the average delta y value of the data to use bins of this width
nrbins = ( max(Z.x_raw(:))-min(Z.x_raw(:)) ) / (dx_avg*2);
Emp = variogram(coord_x,value_x,'nrbins',nrbins,'plotit',false,'maxdist',max(Z.x_raw),'subsample',20000);
assert(all(~isnan(Emp.val)),'problem');

% LSCurve fit
myfun = @(x,h) semivariogram1D(h,c(1),x,'sph', c(2));

range0 = 100;
x_range = lsqcurvefit(myfun,range0,Emp.distance,Emp.val,[],[],options);

% plotit
if plotit
    subplot(1,2,2);hold on;
    plot(Emp.distance,Emp.val,'--o')
    plot(Emp.distance,myfun(x_range,Emp.distance))
    legend('Empirical Variogram vertical','Fitted Variogram horizontal')
end


%% Return

range = [x_range y_range];
disp(['The fitted range is: ', num2str(range)])
disp(['The fitted  is: ', num2str(c)])
keyboard
end



%% Manually Calculated
% Cal.distance = 0:0.1:max(Emp.distance);
% Cal.c0 = 1;
% Cal.model = 'sph';
% Cal.nugget = 0;
% Cal.range = 4.4;
%
% Cal.val = semivariogram1D(Cal.distance,Cal.c0,Cal.range,Cal.model,Cal.nugget);
% figure;hold on;
% plot(Emp.distance,Emp.val,'--o')
% plot(Cal.distance,Cal.val);
% ylabel('\gamma(h)'); xlabel('h')
% legend('Empirical Variogram vertical','Empirical Variogram horizontal')


%% OTHER OLD STUFF
% gy=Emp.val(~isnan(Emp.val(:,1)),1);
% hy=Emp.distance(~isnan(Emp.val(:,1)));
% gx=Emp.val(~isnan(Emp.val(:,2)),2);
% hx=Emp.distance(~isnan(Emp.val(:,2)));
%
% if any(diff(gx)<0)
%     hx=hx(1:find(diff(gx)<0,1));
%     gx=gx(1:find(diff(gx)<0,1));
% end
% if any(diff(gy)<0)
%     hy=hy(1:find(diff(gy)<0,1));
%     gy=gy(1:find(diff(gy)<0,1));
% end
%
% xdata = [hy  zeros(size(hy)); zeros(size(hx)) hx];
% ydata = [gy; gx];
%
% %
%
%
% figure; hold on;
% plot(xdata,ydata,'o')
% plot(xdata,myfun(x0,xdata),'x')
%
% legend('Empirical Variogram vertical','Empirical Variogram horizontal')
% % scatter(X.x,X.y,[],X.d)
