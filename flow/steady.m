function [p,ux,uy] = steady(perm,P0,P1,viscosite,phi,DX,DY)

% [p,ux,uy] = steady(perm,P0,P1,viscosite,phi);

% resolution du probleme des ecoulements
% conditions: monophasique, incompressible, 2D, kxx = kyy, DX=DY
% differences finies, 5 points
% perm: champ de permeabilite
% P0: pression 
% viscosite: viscosite du fluide
% phi: porosite (constante)
%
% p: champ de pression
% ux et uy: champs de vitesses selon X et Y



% matrice des transmissivites et conditions aux limites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ni,Nj] = size(perm);
T = sparse(Ni*Nj,Ni*Nj);
p = zeros(Ni*Nj,1);
indice_connu = 1;
indice_pasconnu = 1;
rapportX = DY/DX;
rapportY = 1/rapportX;


for i = 1:Ni
  for j = 1:Nj

    k = j+(i-1)*Nj;
    terme1 = 0;
    terme2 = 0;
    terme3 = 0;
    terme4 = 0;
    % diagonale
    if (i < Ni)
      terme1 = rapportX*khar(perm(i,j),perm(i+1,j))/viscosite;
    end
    if i > 1
      terme2 = rapportX*khar(perm(i-1,j),perm(i,j))/viscosite;
    end
    if j < Nj
      terme3 = rapportY*khar(perm(i,j),perm(i,j+1))/viscosite;
    end
    if j > 1
      terme4 = rapportY*khar(perm(i,j-1),perm(i,j))/viscosite;
    end
    T(k,k) = -terme1-terme2-terme3-terme4;

    % bande 1
    if i < Ni
      b1 = j+i*Nj;
      T(k,b1) = terme1;
    end    

    % bande 2
    if i > 1
      b2 = j+(i-2)*Nj;
      T(k,b2) = terme2;
    end


    % bande 3
    if j < Nj
      b3 = j+1+(i-1)*Nj;
      T(k,b3) = terme3;
    end


    % bande 4
    if j > 1
      b4 = j-1+(i-1)*Nj;
      T(k,b4) = terme4;
    end

    % Points pour lesquels on connait la pression
    if j == 1
      p(k) = P0;
      connu(indice_connu) = k;
      indice_connu = indice_connu+1;
    elseif j == Nj
      p(k) = P1;
      connu(indice_connu) = k;
      indice_connu = indice_connu+1;
    else
      p(k) = 0;
      pasconnu(indice_pasconnu) = k;
      indice_pasconnu = indice_pasconnu+1;
    end

  end
end


% systeme d'equations a resoudre
% extraction des pressions connues
TT = T(pasconnu,pasconnu);
% conversion en sparse matrix
SPT = sparse(TT);
n = length(pasconnu);
C = zeros(n,1);
j = 1;
for i = 1:n
  C(j) = -T(pasconnu(i),:)*p;
  j = j+1;
end


% pressions
p(pasconnu) = SPT\C;
clear C SPT T TT


% vitesses en x
ux = zeros(Ni,Nj-1);
for j = 1:Nj-1
  for i = 1:Ni
    ux(i,j) = -rapportX*khar(perm(i,j),perm(i,j+1))*(p(j+1+(i-1)*Nj)-p(j+(i-1)*Nj))/(viscosite*phi*DX);
  end
end


% vitesses en y 
uy = zeros(Ni-1,Nj);
for j = 1:Nj
  for i = 1:Ni-1
    uy(i,j) = -rapportY*khar(perm(i,j),perm(i+1,j))*(p(j+i*Nj)-p(j+(i-1)*Nj))/(viscosite*phi*DY);
  end
end


p = reshape(p,Nj,Ni)';