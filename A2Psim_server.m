addpath(genpath('./.')); 
load('result-A2PK/GEN-Run_1_2017-05-07_14-37');

% Nscore of Sigma. 
[kern.prior,kern.axis_prim] = ksdensity(sigma_true(:));
parm.nscore=1;
Nscore = nscore(kern, parm, 0);
NSigma.x=Sigma.x_raw; NSigma.y=Sigma.y_raw; 
[NSigma.X, NSigma.Y] = meshgrid(Sigma.x_raw, Sigma.y_raw);
d = reshape(Nscore.forward(Sigma.d_raw(:)) ,numel(NSigma.y),numel(NSigma.x));
zt = reshape(Nscore.forward(sigma_true(:)), grid_gen.ny, grid_gen.nx);

% Built G
G = zeros(numel(NSigma.X),grid_gen.nxy);
for ij=1:grid_gen.nxy
     [~,id_z]=min((grid_gen.X(ij)-NSigma.X(:)).^2 + (grid_gen.Y(ij)-NSigma.Y(:)).^2);
     G(id_z,ij)=1;%/(NSigma.dx*NSigma.dy);
end
for ij = 1:numel(NSigma.X)
    G(ij,G(ij,:)==1) = 1 ./ sum(G(ij,:)==1);
end


% Compare actual ERT inverted and corresponding averaged Electricit
dG = reshape(G * zt(:), numel(NSigma.y), numel(NSigma.x));
Gt = d(:) * zt(:)' / (zt(:) * zt(:)');
dGt = reshape(Gt * zt(:), numel(NSigma.y), numel(NSigma.x));

% kriging prediction
Cz = covardm_perso([grid_gen.X(:) grid_gen.Y(:)], [grid_gen.X(:) grid_gen.Y(:)], gen.covar);
Czd = Cz * G';
Cd = G * Czd;
W=zeros(grid_gen.nxy,numel(NSigma.X));

parpool(48);
parfor ij=1:grid_gen.nxy
    W(ij,:) = Cd \ Czd(ij,:)';
end

zh = reshape( W * d(:), grid_gen.ny,grid_gen.nx);

save('result-A2PK/A2Psim-Run_1_2017-05-07_14-37_2','Cd','Czd','W','G','Gt','d','zt','zh','dG','dGt');

