function [temps,nombre] = tracertimes(time)

% [temps,nombre] = tracertimes(time);
% courbes des temps d'arrivee

[nn,temps] = hist(time,100);
m = length(nn);
nn(3:m+2) = nn(1:m);
temps(3:m+2) = temps(1:m);
nombre(1) = 0;
temps(1) = 0;
temps(2) = floor(temps(3));
nombre(2) = 0;
nombre(3) = nn(3);
for i = 4:m+2
  nombre(i) = nombre(i-1)+nn(i);
end
nombre = nombre./nombre(length(nn));
clear nn