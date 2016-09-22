function parm = kriginginitiaite(parm)

if isfield(parm, 'gen')
    if isfield(parm.gen.covar, 'modele') && isfield(parm.gen.covar, 'c')
        % old version of gen
        parm.k.covar=struct('model',num2cell(parm.gen.covar.modele(:,1)),...
                'c0',num2cell(parm.gen.covar.c)',...
                'range',num2cell(parm.gen.covar.modele(:,2:3)',1)'...
                ,'azimuth',num2cell(parm.gen.covar.modele(:,4)));
    elseif isfield(parm.gen, 'covar')
        % new version of gen
        assert(isfield(parm.gen.covar, 'model'),'parm.covar is not properly define')
        assert(isfield(parm.gen.covar, 'c0'),'parm.covar is not properly define')
        assert(isfield(parm.gen.covar, 'range'),'parm.covar is not properly define')
        assert(isfield(parm.gen.covar, 'azimuth'),'parm.covar is not properly define')
        parm.k.covar    = parm.gen.covar;
    else
        error('No parm.gen readable')
    end
elseif isfield(parm, 'k') && isfield(parm.k, 'model') && isfield(parm.k, 'c')
    % old version
    parm.k.covar=struct('model',num2cell(parm.k.model(:,1)),...
        'c0',num2cell(parm.k.c)',...
        'range',num2cell(parm.k.modele(:,2:3)',1)'...
        ,'azimuth',num2cell(parm.k.model(:,4)));
else
    error('You need to define parm.covar or parm.gen.covar')
    

end
clear parm.k.model parm.k.c

for i=1:numel(parm.k.covar)
    
    % Patch for old version
    switch parm.k.covar(i).model
        case 1
            parm.k.covar(i).model='nugget';
        case 2
            parm.k.covar(i).model='exponential';
        case 3
            parm.k.covar(i).model='gaussian';
        case 4
            parm.k.covar(i).model='spherical';
    end

    switch parm.k.covar(i).model
        case 'nugget'
            parm.k.covar(i).g = @(h,r) h==0;
            intvario=1;
        case 'triangle'
            assert(numel(parm.gen.covar.range)==1,'only valid in 1D')
            intvario=1;
            parm.k.covar(i).g = @(h) max(1-h,0);
        case 'circular'
            intvario=1.17;
            parm.k.covar(i).g = @(h) 2/pi*(acos(min(h,1))-min(h,1).*sqrt(1-min(h,1).^2));
        case 'spherical'
             intvario=1.3;
            parm.k.covar(i).g = @(h) 1-3/2*min(h,1)+1/2*min(h,1).^3;
        case 'cubic'
            intvario=1.43;
            parm.k.covar(i).g = @(h) 1 - 7*min(h,1).^2 + 35/4*min(h,1).^3 - 7/2*min(h,1).^5 + 3/4*min(h,1).^7;
        case 'exponential'
            intvario = .41;
            parm.k.covar(i).g = @(h) exp(-h);
        case 'gaussian'
            intvario=.58;
            parm.k.covar(i).g = @(h) exp(-h.^2);
        case 'stable'
            assert(isfield(parm.gen.covar, 'alpha'),'alpha covar is not properly define')
            assert(parm.gen.covar(i).alpha>0 && parm.gen.covar(i).alpha<=2,'Alpha value not possible')
            intvario=.41;
            parm.k.covar(i).g = @(h) exp(-(h).^parm.gen.covar(i).alpha);
        case 'power'
            assert(isfield(parm.gen.covar, 'alpha'),'alpha covar is not properly define')
            assert(parm.gen.covar(i).alpha>0 && parm.gen.covar(i).alpha<2,'Alpha value not possible')
            parm.k.covar(i).g = @(h) 1-h.^parm.gen.covar(i).alpha;
            warning('Approx the integrale')
            intvario=1.5;
        case 'k-bessel'
            assert(isfield(parm.gen.covar, 'alpha'),'alpha covar is not properly define')
            assert(parm.gen.covar(i).alpha>0 && parm.gen.covar(i).alpha<2,'Alpha value not possible')
            intvario=[.35 .5];
            parm.k.covar(i).g = @(h) 1/(2^(parm.gen.covar(i).alpha-1) * gamma(parm.gen.covar(i).alpha)) .* max(h,eps).^parm.gen.covar(i).alpha .* besselk(parm.gen.covar(i).alpha,max(h,eps));
        case 'logarithmic'
            parm.k.covar(i).g = @(h) 1-log(h+1);
            warning('Approx the integrale')
            intvario=.7;
        case 'cauchy'
            assert(isfield(parm.gen.covar, 'alpha'),'alpha covar is not properly define')
            assert(parm.gen.covar(i).alpha>0,'Alpha value not possible')
            parm.k.covar(i).g = @(h) (1+h.^2).^parm.gen.covar(i).alpha;
            warning('Approx the integrale')
            intvario=1;
        case 'hyperbolic'
            parm.k.covar(i).g = @(h) 1./(1+h);
            warning('Approx the integrale')
            intvario=[.2 .05];
        case 'cardinal sine'
            intvario=.2;
            parm.k.covar(i).g = @(h) sin(max(eps,h))./max(eps,h);
        case 'matheron'
            parm.k.covar(i).g = @(h) 1./(h);
        otherwise
            error('Variogram type not defined')
    end
    
    if numel(parm.gen.covar(1).range)==1 || numel(intvario)==1
        parm.k.covar(i).range=parm.gen.covar(i).range*intvario(1);
    elseif numel(parm.gen.covar(1).range)==2
        parm.k.covar(i).range=parm.gen.covar(i).range*intvario(2);
    elseif numel(parm.k.covar(i).azimuth)==3
        parm.k.covar(i).range=parm.gen.covar(i).range*intvario(3);
    end
    
end

for i=1:numel(parm.k.covar)
    if numel(parm.gen.covar(1).range)==1 || numel(parm.gen.covar(1).azimuth)==0
        parm.k.covar(i).cx = 1/diag(parm.k.covar(i).range(1));
    elseif numel(parm.gen.covar(1).range)==2
        ang=parm.k.covar(i).azimuth; cang=cos(ang/180*pi); sang=sin(ang/180*pi);
        rot = [cang,-sang;sang,cang];
        parm.k.covar(i).cx = rot/diag(parm.k.covar(i).range);
    elseif numel(parm.k.covar(i).azimuth)==3
        error('3D')
    end
end

end