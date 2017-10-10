function [x fval history] = myproblem(x0)



history = {};
options = optimset('MaxFunEvals',3*24,'OutputFcn',@myoutput,'Display','iter');
[x fval] = fminsearch(@objfun, x0,options);


    function stop = myoutput(x,optimValues,state)
        optimValues.state = state;
        optimValues.x = x;
        stop = false;
        if isequal(state,'iter')
            history = {history{:}, optimValues};
            save(['result-BSS/Cal00_', num2str(optimValues.iteration) ,'.mat'],'history');
        end
        switch state
            case 'init'
                %header = ' Iteration   Func-count     min f(x)      x        Procedure';
                %disp(header)
                disp('Starting')
            case 'iter'
                disp(['---------- iter: ', num2str(optimValues.iteration), ' ------'])
                disp(optimValues.x)
                disp(optimValues.fval)
                %fprintf(' %5.0f        %5.0f     %12.6g    %12.6g    %s\n', optimValues.funccount, optimValues.fval, x, optimValues.procedure);
            case 'done'
                disp('Done');
                
        end
    end




    function z = objfun(x)
        
        load('result-BSS/GEN-Run_1_2017-05-07_14-37','gen','grid_gen','K','Sigma','kern');
        [Sigma.X,Sigma.Y]=meshgrid(Sigma.x, Sigma.y);
        addpath(genpath('./.'));
        parm.k.covar = gen.covar;

        parm.unit='';
        parm.nscore = 1;
        parm.par = 1;
        parm.par_n = 48;

        parm.k.nb = [0 0 0 0 0; 10 10 10 10 20];
        parm.cstk = 1;
        parm.seed = 'shuffle';
        parm.scale=[grid_gen.sx;grid_gen.sy]; % no multigrid
        parm.saveit = false;
        parm.k.method = 'sbss'; % sort, sbss (29.6) or smart minkmex
        parm.k.quad = 0;
        parm.k.wradius = 1.3;
        parm.plot.kernel=0;
        parm.plot.ns= 1;
        parm.plot.krig=0;
        parm.path='quasirandom';
        parm.path_random=1;
        parm.kernel=kern;


        Klog=K;
        Klog.d=log(Klog.d);

        parm.n_realisation  = parm.par_n*4;
        parm.aggr.method='rad';
        parm.aggr.A=x(1);
        parm.aggr.B=x(2);

        [Res,~,kern,~,~,~] =  BSGS(Klog,Sigma,grid_gen,parm);

        Res_m_ns = Res.m_ns;
        id = Res.x<parm.k.covar(1).range0(1).*parm.k.wradius;
        Gamma_t = (1-parm.k.covar(1).g(grid_gen.x/parm.k.covar(1).range(1)))';
        Gamma_t_id = Gamma_t( any( bsxfun(@eq, grid_gen.x,Res.x(id)) ));
        XY = kern.XY;
        Sigma_d = Sigma.d(:);
        Res_m = Res.m;
        dens = parm.kernel.dens(:);


        E1 = nan(sum(id),parm.n_realisation);
        E2 = nan(numel(parm.kernel.dens),parm.n_realisation);
        parfor i_realisation=1:parm.n_realisation
            gamma_x = variogram_gridded_perso(Res_m_ns{i_realisation});
           E1(:,i_realisation) = gamma_x(id)-Gamma_t_id;
           E2(:,i_realisation) = ksdensity([Sigma_d Res_m{i_realisation}(:)],XY)-dens;
        end
        OF1 = sqrt(mean(mean(E1,2).^2));
        OF2 = sqrt(mean(mean(E2,2).^2));

        z = OF1/.35 + OF2/.015;

    end
end
