function res = OF_fx(Klog,Sigma,grid_gen,parm,id,Gamma_t,XY,jpdf,Sigma_d,a,b)

parm.aggr.A=a;
parm.aggr.B=b;


Res =  BSGS(Klog,Sigma,grid_gen,parm);

[OF1, ~] = OF_fx_2(Res,id,Gamma_t,XY,jpdf,Sigma_d,parm);

res = [OF1];%, OF2];

end

function [OF1, OF2] = OF_fx_2(Res,id,Gamma_t,XY,jpdf,Sigma_d,parm)


E1 = nan(sum(id),parm.n_realisation);
%E2 = nan(numel(jpdf),parm.n_realisation);
parfor i_realisation=1:parm.n_realisation
    gamma_x = variogram_gridded_perso(Res.m_ns{i_realisation});
    E1(:,i_realisation) = gamma_x(id);
    %E2(:,i_realisation) = ksdensity([Sigma_d Res.m{i_realisation}(:)],XY)-jpdf;
end
-Gamma_t
OF1= sqrt(mean(mean(E1,2).^2));
OF1= (numel(Gamma_t)-1) * (Gamma-Gamma_t) / S_Gamma * (Gamma-Gamma_t)';
OF2=NaN;% sqrt(mean(mean(E2,2).^2));
end




