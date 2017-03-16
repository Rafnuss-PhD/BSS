function OF = OF(filename)

load(filename)




% 1st 2nd order stat.
OF_true.mean = mean(Prim_true(:));
OF_true.std = std(Prim_true(:));
% Variogram
[OF_true.vario_x_field, OF_true.vario_y_field] = variogram_gridded_perso(Prim_true);
OF_true.vario_x_theo = semivariogram1D(grid{end}.x,1,parm.gen.covar.modele(1,2),parm.gen.covar.modele(1),0);
OF_true.vario_y_theo = semivariogram1D(grid{end}.y,1,parm.gen.covar.modele(1,3),parm.gen.covar.modele(1),0);

OF.hist_pts = -4:.1:4;

% Entropy
OF.H_true = 1/2*log(2*pi)+1/2; % crossentropy(normrnd(0,1,100000,1));


% Parameter for the ERT inversion
dmin = parm.gen.Rho.dmin;
grid_Rho = parm.gen.Rho.grid;
elec = parm.gen.Rho.elec;
dmin.filepath       = 'result-time-log-calibration/4-Calibration/IO-file/';
dmin.header         = 'Forward';  % title of up to 80 characters
dmin.job_type       = 0;
dmin.readonly       = 0;



for i_realisation=1:parm.n_realisation
    
    %weight
    for i_scale=1:numel(Res)
        OF.w(i_realisation,i_scale) = Res{i_scale}.w_krig;
    end
    
    % 1st 2nd order stat.
    OF.mean(i_realisation) = mean(Res{end}.m_ns{i_realisation}(:))';
    OF.std(i_realisation) = std(Res{end}.m_ns{i_realisation}(:))';
    OF.hist(i_realisation,:) = ksdensity(Res{end}.m_ns{i_realisation}(:),OF.hist_pts);
    
    % Kullback–Leibler divergence
    OF.H(i_realisation) = crossentropy(normpdf(OF.hist_pts,0,1),Res{end}.m_ns{i_realisation}(:)) - OF.H_true ;
    
    % Varigram
    [OF.gamma_x(:,i_realisation), OF.gamma_y(:,i_realisation)] = variogram_gridded_perso(Res{end}.m_ns{i_realisation});
    
    % ERT
    delete([dmin.filepath '*']);
    rho_true            =  1000./Res{end}.m{i_realisation}; %sigma_true;%
    dmin.rho_min        = min(rho_true(:));
    dmin.rho_avg        = mean(rho_true(:));
    dmin.rho_max        = max(rho_true(:));
    
    f                   = griddedInterpolant({grid{end}.y,grid{end}.x},rho_true,'nearest','nearest');
    dmin.rho_true       = f({grid_Rho.y,grid_Rho.x});
    dmin.filename       = 'gtrue.dat';
    dmin.num_regions    = 1+numel(dmin.rho_true);
    
    dout                = Matlat2R2(grid_Rho,dmin,elec); % write file and run forward modeling
    OF.ert(:,i_realisation) = dout.output.pseudo;
    OF.ert_interp(:,:,i_realisation) = dout.output.pseudo_interp;
    OF.ert_err(i_realisation) = 1/parm.gen.Rho.f.num_ind_meas*sum( (OF.ert(:,i_realisation) - parm.gen.Rho.f.output.pseudo).^2);
    
    
    % True field not knwon
    OF.true_err(i_realisation) = 1/Res{end}.nx/Res{end}.ny * sum(( Res{end}.m{i_realisation}(:)-Prim_true(:) ).^2);
    
end

if 0
    % Figure
    figure(1);
    subplot(3,2,[1 2]); hold on;
    ksdensity(Prim_true(:))
    ksdensity(Res{end}.m{i_realisation}(:))
    axis tight; legend('true','simulation 1')
    
    subplot(3,2,3); hold on;
    plot(grid{end}.x, OF_true.vario_x_theo)
    plot(grid{end}.x, OF_true.vario_x_field)
    plot(grid{end}.x, OF.gamma_x{i_realisation},'.-')
    axis tight; xlim([0 parm.gen.covar.modele(1,2)*1.3])
    
    
    subplot(3,2,4); hold on;
    plot(grid{end}.y, OF_true.vario_y_theo)
    plot(grid{end}.y, OF_true.vario_y_field)
    plot(grid{end}.y, OF.gamma_y{i_realisation},'.-')
    axis tight; xlim([0 parm.gen.covar.modele(1,3)*1.3])
    
    subplot(3,2,[5 6]);
    imagesc(flipud(OF.ert_interp{i_realisation}))
    axis tight
end
end