function [k] = khar(k1,k2);

if k1 ~= 0 & k2 ~= 0
  k = 2/(1/k1+1/k2);
elseif k1 == 0 & k2~= 0
  k = k2;
elseif k1 ~= 0 & k2 == 0
  k = k1;
else
  k = 0;
end