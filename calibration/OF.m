
% 
function OF(filename)

load(filename)

% True Field calculation
% Entropy
OF_true.H = crossentropy(Prim_true(:));
% Variogram
[OF_true.vario_x_field, OF_true.vario_y_field] = variogram_gridded_perso(Prim_true);

OF_true.vario_x_theo = semivariogram1D(grid_gen.x,1,parm.gen.covar.modele(1,2),parm.gen.covar.modele(1),0);
OF_true.vario_y_theo = semivariogram1D(grid_gen.y,1,parm.gen.covar.modele(1,3),parm.gen.covar.modele(1),0);



for i_realisation=1:parm.n_realisation
    
    M = Res{end}.m{i_realisation};
    
    % Jensen-Shannon divergence
    M_H = crossentropy(Prim_true(:),M(:));
    H = M_H - OF_true.H;

    % Varigram
    [gamma_x, gamma_y] = variogram_gridded_perso(M);

    % Field
    error = M(:)-Prim_true(:);
         

    
end

% Figure
figure(1); 
subplot(3,2,[1 2]); hold on;
ksdensity(Prim_true(:))
ksdensity(M(:))
axis tight;

subplot(3,2,3); hold on;
plot(grid_gen.x, OF_true.vario_x_theo)
plot(grid_gen.x, OF_true.vario_x_field)
plot(grid{end}.x, gamma_x,'.-')
axis tight; xlim([0 parm.gen.covar.modele(1,2)*1.3])


subplot(3,2,4); hold on;
plot(grid_gen.y, OF_true.vario_y_theo)
plot(grid_gen.y, OF_true.vario_y_field)
plot(grid{end}.y, gamma_y,'.-')
axis tight; xlim([0 parm.gen.covar.modele(1,3)*1.3])

subplot(3,2,[5 6]); hold on;
plot(Prim_true(:),error,'.')

end