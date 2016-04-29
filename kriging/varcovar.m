function [CY,vario,l,id] = varcovar(Res, CY,vario,ord)
l=Res{end}.lambda;
ny=Res{end}.ny;
nx=Res{end}.nx;
nxy=nx*ny;

if CY
    CY =  (l) \ transpose(inv(l));
else
    CY= [];
end

if vario
    vario{3} = nan(ny*2-1,nx*2-1,nxy);
    for i=1:nxy
        [id_y,id_x]=ind2sub([ny,nx],i);
        vario{3}( (ny-id_y)+(1:ny) , (nx-id_x)+(1:nx),i) = reshape(full(CY(i,:)),ny,nx);
    end
    
    vario{2} = nanmean(vario{3},3);
    vario{1}.h = [ Res{end}.varcovar.vario_2d(Res{end}.ny,Res{end}.nx), mean( [Res{end}.varcovar.vario_2d(Res{end}.ny,Res{end}.nx-1:-1:1); Res{end}.varcovar.vario_2d(Res{end}.ny,Res{end}.nx+1:end)])];
    vario{1}.v = [ Res{end}.varcovar.vario_2d(Res{end}.ny,Res{end}.nx), mean( [Res{end}.varcovar.vario_2d(Res{end}.ny-1:-1:1,Res{end}.nx)'; Res{end}.varcovar.vario_2d(Res{end}.ny+1:end,Res{end}.nx)'])];
else
    vario=[];
end


if ord
    id=[];
    for i_scale=1:parm.n_scale
        id = [id; Res{i_scale}.varcovar_id(Res{i_scale}.sim.xy_r)'];
    end
    
    l = l(id,id);
    if CY
        CY = CY(id,id);
    end
    if vario;
        error('not implemented yet')
    end
end


    