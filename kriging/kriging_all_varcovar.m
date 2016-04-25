function varcovar_lambda=kriging_all_varcovar(Res,parm)


varcovar_lambda = nan(Res.nx*Res.ny,Res.nx*Res.ny);

C=covardm_perso([Res.X(:) Res.Y(:)],[Res.X(:) Res.Y(:)],parm.k.model,parm.k.var,parm.k.cx);


parfor i = 1:Res.nx*Res.ny
    a0_C=[C(1:i-1,i) ; C(i+1:end,i)];
    
    ab_C= [C(1:i-1,1:i-1)  C(1:i-1,i+1:end);...
        C(i+1:end,1:i-1) C(i+1:end,i+1:end)];
    
    pt_krig_lambda = ab_C \ a0_C; % Ordinary
    pt_krig_s = sum(parm.k.var) - pt_krig_lambda'*a0_C;


    varcovar_lambda(i,:) = 1./sqrt(pt_krig_s)  * [-pt_krig_lambda(1:i-1); 1; -pt_krig_lambda(i:end)];
    
    if i/100==1
        disp(i)
    end
end

end
