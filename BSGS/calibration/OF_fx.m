function res = OF_fx(Klog,Sigma,grid_gen,parm,a,b)

parm.aggr.A=a;
parm.aggr.B=b;

[Res,~,kern,~,~,~] =  BSGS(Klog,Sigma,grid_gen,parm);

Res_m_ns = Res.m_ns;
id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));
XY = kern.XY;
Sigma_d = Sigma.d(:);
Res_m = Res.m;
dens = parm.kernel.dens(:);


E1 = nan(sum(id),parm.n_realisation);
E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
parfor i_realisation=1:parm.n_realisation
    gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
   E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
   E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
end
OF1 = sqrt(mean(mean(E1,2).^2));
OF2 = sqrt(mean(mean(E2,2).^2));

res = OF1/.35 + OF2/.015;

end

function [OF1, OF2] = OF_fx_2(Res_m_ns,Res_m,id,Gamma_t_id,XY,Sigma_d,dens,parm)


E1 = nan(sum(id),parm.n_realisation);
E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
parfor i_realisation=1:parm.n_realisation
    gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
   E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
   E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
end
OF1 = sqrt(mean(mean(E1,2).^2));
OF2 = sqrt(mean(mean(E2,2).^2));

end




