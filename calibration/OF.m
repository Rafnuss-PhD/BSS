function res = OF(Klog,Sigma,grid_gen,parm,id,Gamma_t,X,Y,jpdf,Sigmad3D, x)

parm.aggr.fx = @(parm,grid,i_scale,i_pt) x(i_scale);
Res =  BSGS(Klog,Sigma,grid_gen,parm);

[E1, E2] = OF2(Res,id,Gamma_t,X,Y,jpdf,Sigmad3D,grid_gen,parm);

res = [E1 E2];

end

function [E1, E2] = OF2(Res,id,Gamma_t,X,Y,jpdf,Sigmad3D,grid_gen,parm)
% Variogram
for i_realisation=1:parm.n_realisation
    [gamma_x(:,i_realisation), ~] = variogram_gridded_perso(Res{end}.m_ns{i_realisation});
end


E1 = sqrt(mean((mean(gamma_x(id,:),2) - Gamma_t(:)).^2));

% Joint PDF
Res3D = reshape([Res{end}.m{:}], [grid_gen.ny, grid_gen.nx, parm.n_realisation]);
Resjpdf=reshape(ksdensity([Sigmad3D(:) Res3D(:)],[X(:),Y(:)]),size(X,1),size(X,2));
E2 = sqrt(mean((Resjpdf(:) - jpdf(:)).^2));
end




