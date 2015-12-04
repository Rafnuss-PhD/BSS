function [k]=covardm_perso(x,x0,model,c,cx)
% [k]=Covar(x,x0,model,c)
% Fonction pour calculer les covariances avec des mod`eles sp�cifi�es comme avec cokri
% La fonction calcule pour toutes les positions de x (n1) avec toutes les positions de x0 (n2)  K est donc n1xn2
% auteur D. Marcotte avril 2002

% here we define the equations for the various covariograms. Any new model
% can be added here.
% k=[];
% Gam=['h==0                                              '; %nugget
%      'exp(-h)                                           '; %exponential
%      'exp(-(h).^2)                                      '; %gaussian
%      '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)             '; %spherical
%      '1-h                                               '; %linear
%      '1-3*min(h,1).^2+2*min(h,1).^3                     '; %modele Trochu
%      '(h.^2).*log(max(h,eps))                           '; %spline plaque mince
%      '(h.^2+1).^(-0.5)                                  '; %mod�le gravim�trique (Cauchy avec b=0.5)
%      '(h.^2+1).^(-1.5)                                  '; %modele magn�tique (Cauchy avec b=1.5)
%      'sin(max(eps,h*2*pi))./max(eps,h*2*pi)             '; %effet de trou sinusoidal
%      'cos(h*2*pi)                                       '; %effet de trou cosinusoidal
%      '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)+1-h         '];%spherique+lineaire


% definition of some constants


%  n1=size(x,1); % d dimension de l'espace
%  n2=size(x0,1);
% r=size(c,1);
% assume d=2, p=1;
% r=rp/p;  % nombre de structures
% cx=[x(:,1:d);x0];
% nm=size(model,2);

% ne pas permettre des port�es de 0 en input
% assert(all(model(:,2)>0),'range needs to be greater than 0')

% calculer les covariances

k=zeros(size(x,1),size(x0,1));

for i=1:numel(c)
    % calculation of matrix of reduced rotated distances H
    % [t1]=trans(x,model,i);
    % [t2]=trans(x0,model,i);
    %     t1 = x*cx{i};
    %     t2 = x0*cx{i};
    
    %    h=0;
    %    for id=1:d
    %       h=h+(t1(:,id)*ones(1,n2)-ones(n1,1)*t2(:,id)').^2;
    %    end
    %    h=sqrt(h);
    if size(x0,1)==1
        h=sqrt(sum(bsxfun(@minus,x*cx{i},x0*cx{i}).^2,2));
    elseif all(size(x)==size(x0)) && all(x(:)==x0(:))
        h=squareform(pdist(x*cx{i}));
    else
        error('should be here...')
        h=pdist2(x*cx{i}, x0*cx{i});
    end
    
    
    
    % evaluation of the current basic structure
    switch model(i,1)
        case 1
            g = h==0; %nugget
        case 2
            g = exp(-h)                                           ; %exponential
        case 3
            g = exp(-(h).^2)                                      ; %gaussian
        case 4
            g = 1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)             ; %spherical
        case 5
            g = 1-h                                               ; %linear
        case 6
            g = 1-3*min(h,1).^2+2*min(h,1).^3                     ; %modele Trochu
        case 7
            g = (h.^2).*log(max(h,eps))                           ;%spline plaque mince
        case 8
            g = (h.^2+1).^(-0.5)                                  ; %mod�le gravim�trique (Cauchy avec b=0.5)
        case 9
            g = (h.^2+1).^(-1.5)                                  ; %modele magn�tique (Cauchy avec b=1.5)
        case 10
            g = sin(max(eps,h*2*pi))./max(eps,h*2*pi)             ; %effet de trou sinusoidal
        case 11
            g = cos(h*2*pi)                                       ; %effet de trou cosinusoidal
        case 12
            g = 1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)+1-h;         %spherique+lineaire
    end
    % g=eval(Gam(model(i,1),:));
    %ji=(i-1)*2; js=i;
    %k=k+kron(g,c(ji:js,:));
    
    k=k+kron(g,c(i));
end

% function [cx,rot]=trans(cx,model,im)
% function [cx,rot]=trans(cx,model,im);
%
% TRANS is called from COKRI2. It takes as input original coordinates and
%       return the rotated and reduced coordinates following specifications
%       described in model(im,:)
% Rotations are all performed anticlockwise with the observer located on the positive side of
% the axis and looking toward the origin. In 3D, rotations are performed first along z,
% then along rotated y and then along twice rotated x.
% Author: D. Marcotte
% Version 2.1  97/aug/18

% some constants are defined

% assume d)=2, p=4
% [~,d]=size(cx);
% [~,p]=size(model);

% check for 1-D or isotropic model

% if p-1>d,

% perform rotation counterclockwise

%   if d==2,
%       ang=model(im,4); cang=cos(ang/180*pi); sang=sin(ang/180*pi);
%       rot=[cang,-sang;sang,cang];
%    else
%
%       % rotation matrix in 3-D is computed around z, y and x in that order
%
%       angz=model(im,7); cangz=cos(angz/180*pi); sangz=sin(angz/180*pi);
%       angy=model(im,6); cangy=cos(angy/180*pi); sangy=sin(angy/180*pi);
%       angx=model(im,5); cangx=cos(angx/180*pi); sangx=sin(angx/180*pi);
%       rotz=[cangz,-sangz,0;sangz,cangz,0;0 0 1];
%       roty=[cangy,0,sangy;0 1 0;-sangy,0,cangy];
%       rotx=[1 0 0;0 cangx -sangx;0 sangx cangx];
%       rot=rotz*roty*rotx;
%    end

% rotation is performed around z, y and x in that order, the other coordinates are left unchanged.

% dm=min(3,d);
%cx(:,1:dm)=cx(:,1:dm)*rot;
% cx=cx*rot;
% t=[model(im,2:1+dm),ones(d-dm,1)];
% t=diag(t);
% else
%    t=eye(d)*model(im,2);
% end

% perform contractions or dilatations (reduced h)

% cx=cx/t;
%  cx = cx*rot/diag(model(im,2:3));

