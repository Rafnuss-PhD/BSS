function [out, rout] = fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_t_id,Sigma_d,XY,dens,hd)

Rest = nan(ny,nx,parm.n_real);

parm.aggr.T = T; % w_min w_max, T_1, T_2
parm_seed_U = 'shuffle';

parfor i_real=1:parm.n_real
    Res=nan(ny,nx);
    Res(hd.id) = hd.d;
    rng(parm_seed_U);
    U=rand(ny,nx);
    for i_scale = 1:sn
        for i_pt = start(i_scale)+(1:nb(i_scale))
            n = ~isnan(NEIGH(i_pt,:));
            fkrig = normpdf(sec.axis, LAMBDA(i_pt,n)*Res(NEIGH(i_pt,n))', sqrt(S(i_pt)));
            fkrig = fkrig' ./ sum(fkrig);
            
            fsec = sec.pdf(:,path(i_pt));
            
            w = aggr_fx(i_real,i_pt/sum(nb),parm.aggr,i_scale/sn);
            
            
            fa = f0.^0 .* fkrig.^(1-w) .* fsec.^w;
            fa=fa./sum(fa);
            
            cfa = cumsum([0 ; fa(2:end-1)+eps ; eps]) ./ (sum(fa(2:end-1)+eps));
            Res(path(i_pt)) =  interp1(cfa, sec.axis, U(i_pt),'pchip');
            
        end
    end
    Rest(:,:,i_real) = Res;
end


E1 = nan(sum(id),parm.n_real);
E2 = nan(numel(kern.dens),parm.n_real);
parfor i_real=1:parm.n_real
   r = Rest(:,:,i_real);
   gamma_x = variogram_gridded_perso(r);
   E1(:,i_real) = gamma_x(id)-Gamma_t_id;
   dens_r =ksdensity([Sigma_d r(:)],XY);
   E2(:,i_real) = dens_r./sum(dens_r) - dens;
end


OF1=nan(size(parm.aggr.T,1),1);
OF2=nan(size(parm.aggr.T,1),1);
OF1r=nan(size(parm.aggr.T,1),2);
OF2r=nan(size(parm.aggr.T,1),2);
parfor i_t = 1:size(parm.aggr.T,1)
    ii_t=i_t:size(parm.aggr.T,1):parm.n_real;
    OF1(i_t) = sqrt(mean(mean(E1(:,ii_t),2).^2));
    OF2(i_t) = sqrt(mean(mean(E2(:,ii_t),2).^2));
    ii_t2 = round(numel(ii_t)/2);
    OF1r(i_t,:) = [sqrt(mean(mean(E1(:,ii_t(1:ii_t2)),2).^2)) sqrt(mean(mean(E1(:,ii_t(ii_t2+1:end)),2).^2))];
    OF2r(i_t,:) = [sqrt(mean(mean(E2(:,ii_t(1:ii_t2)),2).^2)) sqrt(mean(mean(E2(:,ii_t(ii_t2+1:end)),2).^2))];
end

OF1_range=[.01 .36];
OF2_range=[3.4 10]*10^(-5);

out1 = (OF1-OF1_range(1))./(range(OF1_range));
out2 = (OF2-OF2_range(1))./(range(OF2_range));

out = out1+out2;
rout = range(OF1r)/range(OF1_range) + range(OF2r)/range(OF2_range);

% out = exp( -3*() );


end


function w = aggr_fx(i_real, x, aggr, s)


i_t = mod(i_real,size(aggr.T,1));
if i_t==0;
    i_t=size(aggr.T,1);
end
assert(x<=1,'error')

switch aggr.method
    case 'cst'
        w=aggr.T(i_t);
    case 'step'
        if (x<aggr.T(i_t))
            w=0;
        else 
            w=1;
        end
    case 'linear'
        if (x<aggr.T(i_t,1))
            w  = aggr.T(i_t,3);
        elseif (x<aggr.T(i_t,2))
            w =  aggr.T(i_t,3) + ( x - aggr.T(i_t,1) )/(aggr.T(i_t,2)-aggr.T(i_t,1)) * (aggr.T(i_t,4)-aggr.T(i_t,3));
        else 
            w = aggr.T(i_t,4);
        end
    case 'sigmoid'
        a = aggr.T(i_t,1);
        b = aggr.T(i_t,2);
        w = (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
    case 'mg'
        i_s = ceil(s*size(aggr.T,2));
        w = aggr.T(i_t,i_s);
    otherwise
        error('no Aggr method defined')  
end
end
