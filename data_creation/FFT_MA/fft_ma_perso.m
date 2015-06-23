function A=fft_ma_perso(grid,moy,var,sill,range_x,range_y,rotation)

if ~exist('fft_ma.m','file')
    run /home/raphael/Dropbox/PhD/Bibliography/Software/mGstat/mgstat_set_path.m
end

options.gmean=moy; 
options.gvar=var;
V.par1 = sill;
V.par2 = [range_x rotation range_x/range_y];
V.type='Sph';
V.itype=1;

A=fft_ma(grid.x,grid.y,format_variogram(V),options);

end