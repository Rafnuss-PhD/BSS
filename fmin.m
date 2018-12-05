function [Rest, out, OF1, OF2, E1_x, E1_y, E2] = fmin(T,parm,ny,nx,sn,start,nb,LAMBDA,NEIGH,S,sec,path,f0,id,kern,Gamma_id,Sigma_d,XY,dens,hd)
%[OF1,OF2,OF1std,OF2std]
%[Rest,E1_x, E1_y, E2]
%[out,rout]

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
            
            
            fa = f0.^(1-w(1)-w(2)) .* fkrig.^w(1) .* fsec.^w(2);
            fa=fa./sum(fa);
            
            cfa = cumsum([0 ; fa(2:end-1)+eps ; eps]) ./ (sum(fa(2:end-1)+eps));
            Res(path(i_pt)) =  interp1(cfa, sec.axis, U(i_pt),'pchip');
            
        end
    end
    Rest(:,:,i_real) = Res;
end


E1_x = nan(numel(id.x),parm.n_real);
E1_y = nan(numel(id.y),parm.n_real);
E2 = nan(numel(kern.dens),parm.n_real);

parfor i_real=1:parm.n_real
   r = Rest(:,:,i_real);
   % r = (r-mean(r(:)))/std(r(:));
   [gamma_x, gamma_y] = variogram_gridded_perso(r);
   E1_x(:,i_real) = gamma_x(id.x)-Gamma_id.x(:);
   E1_y(:,i_real) = gamma_y(id.y)-Gamma_id.y(:);
   dens_r =ksdensity([Sigma_d r(:)],XY);
   E2(:,i_real) = dens_r./sum(dens_r) - dens;
end


OF1=nan(size(parm.aggr.T,1),1);
OF2=nan(size(parm.aggr.T,1),1);
%  OF1std=nan(size(parm.aggr.T,1),1);
%  OF2std=nan(size(parm.aggr.T,1),1);
parfor i_t = 1:size(parm.aggr.T,1)
    ii_t=i_t:size(parm.aggr.T,1):parm.n_real;
    OF1(i_t) = sqrt(sum((((1-Gamma_id.x(:)).*mean(E1_x(:,ii_t),2)).^2))) + sqrt(sum((((1-Gamma_id.y(:)).*mean(E1_y(:,ii_t),2)).^2))) ;
    OF2(i_t) = sqrt(mean(mean(E2(:,ii_t),2).^2));
%     ii_t2 = round(numel(ii_t)/2);
%     r1 = [sqrt(mean(mean(E1(:,ii_t(1:ii_t2)),2).^2)) sqrt(mean(mean(E1(:,ii_t(ii_t2+1:end)),2).^2))];
%     OF1std(i_t) = sqrt(sum((r1-OF1(i_t)).^2))
%     r2 = [sqrt(mean(mean(E2(:,ii_t(1:ii_t2)),2).^2)) sqrt(mean(mean(E2(:,ii_t(ii_t2+1:end)),2).^2))];
%     OF2std(i_t) = sqrt(sum((r2-OF2(i_t)).^2))
end

OF1_range= [0.0141 1.6420];
OF2_range=1.0e-04 *[0.3398 0.8857];

out1 = (OF1-OF1_range(1))./(range(OF1_range));
out2 = (OF2-OF2_range(1))./(range(OF2_range));


out = out1+out2;
%rout = OF1std./range(OF1_range) + OF2std./range(OF2_range);

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
        w=aggr.T(i_t,:);
    case 'step'
        if (x<aggr.T(i_t,1))
            w=[0 aggr.T(i_t,2)];
        else 
            w=[aggr.T(i_t,2) 0];
        end
    case 'linear'
        if x<.25
            w  = aggr.T(i_t,1);
        elseif x<.5
            w  = aggr.T(i_t,2);
        elseif x<.75
            w  = aggr.T(i_t,3);
        else
            w  = aggr.T(i_t,4);
        %elseif (x<aggr.T(i_t,2))
        %    w =  aggr.T(i_t,3) + ( x - aggr.T(i_t,1) )/(aggr.T(i_t,2)-aggr.T(i_t,1)) * (aggr.T(i_t,4)-aggr.T(i_t,3));
%         else 
%             w = 0;%aggr.T(i_t,4);
        end
    case 'sigmoid'
        a = aggr.T(i_t,1);
        b = aggr.T(i_t,2);
        w = (atan(a*b) - atan(b*(a -  x )))/(atan(a*b) - atan(b*(a - 1)));
    case 'mg'
        i_s = ceil(s*size(aggr.T,2));
        w = aggr.T(i_t,i_s);
    case 'hybrid'
        range = sum(x>aggr.step)+1;
        w=[aggr.T(i_t,end)-aggr.T(i_t,range) aggr.T(i_t,range)];
    otherwise
        error('no Aggr method defined')  
end
end
