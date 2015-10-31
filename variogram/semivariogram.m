function gamma=semivariogram(h,c0,xrange,yrange,model,n)

H=sqrt( (h(:,1)/xrange).^2 + (h(:,2)/yrange).^2 );

switch model
    case 'sph'
        gamma = c0*ones(size(H));
        gamma(H<1) = c0*( 3/2*H(H<1) - 1/2*H(H<1).^3);
    case 'exp'
        gamma = c0 *( 1 - exp(-3*H));
    case 'gaus'
        gamma = c0 *( 1 - exp(-3*H.^2));
end

% nugget effect
gamma = n+gamma(H~=0);
end
