function [f,k,temps,nombre] = mygradual(gradualparam,Z1,Z2,model,proportions,faciesvalues,xta,yta,xx,moyenne,variance,portee,type,coefnugget,xdata,ydata,welldata,NX,DX,NY,DY,P0,P1,viscosite,phi,Nparticles,datatime,datanombre);


% gradual deformation of the Gaussian white noises
Z = Z1.*cos(pi*gradualparam)+Z2.*sin(pi*gradualparam);

% Gaussian realization constrained to the static data
[k] = domodel(Z,model,proportions,faciesvalues,xta,yta,xx,moyenne,variance,portee,type,coefnugget,xdata,ydata,welldata,NX,DX,NY,DY);

% fluid flow simulation
[p,ux,uy] = steady(k,P0,P1,viscosite,phi,DX,DY);
clear p

% tracer experiment
if P0 > P1
    xo = DX;
else
    xo = (NX-1)*DX;
end
for i = 1:Nparticles
    stop = 0;
    x = [];
    while length(x) == 0 & stop < 10
        yo = rand(1)*NY*DY;
        [x,y,tof] = addstreamlineFD_Flowv2(xo,yo,NX,NY,DX,DY,ux,uy,P0,P1);
        stop = stop+1;
    end
    clear x y
    time(i) = tof;
end
[temps,nombre] = tracertimes(time);
clear time

% objective function
interpnombre = interp1(temps,nombre,datatime);
i = isnan(interpnombre);
indexzero = find(i == 0);
if length(indexzero) > 0
    minindex = indexzero(1)-1;
    n = length(indexzero);
    maxindex = n+1;
    
    indexsup = i;
    if minindex > 0
        indexsup(1:minindex) = 0;
    end
    [index] = find(indexsup == 1);
    interpnombre(index) = 1;
    clear indexsup 
    
    indexinf = i;
    n = length(i);
    indexinf(maxindex:n) = 0;
    [index] = find(indexinf == 1);
    interpnombre(index) = 0;
    clear index indexinf
end
f = .5*sum((interpnombre-datanombre).*(interpnombre-datanombre));

if isnan(f) == 1
    interpnombre
    f
    pause
end
    
clear interpnombre
