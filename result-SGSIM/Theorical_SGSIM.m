clear all; clc;

% Declare variable

n=5;
neigh = 2;
g = [1; sym('g%d', [n 1])];
syms h rr
f = {
    symfun(exp(-h/rr)                                      ,[h,rr]), 2, 'exponential';
    symfun(exp(-(h/rr).^2)                                 ,[h,rr]), pi^.5, 'gaussian';
    symfun(1 - 1.5*h/rr + .5*(h/rr).^3                     ,[h,rr]), 3/5, 'spherical';
    symfun(exp(-sqrt(h/rr))                                ,[h,rr]), 2*gamma((.5+1)/.5), 'stable';
    symfun(1./(1+h/rr)                                     ,[h,rr]), 1/4, 'hyperbolic';
    symfun(sin(h/rr)./(h/rr)                               ,[h,rr]) pi, 'cardinal sine';
    };

% Row-by-row
gammat = [g(end:-1:1); g(2:end)];
for u=1:numel(neigh)
    L{u} = sym('L%d%d', [n n]);
    for i=1:n
        nn=min(i-1,neigh(u));
        Cij = sym('C%d%d', [nn nn]);
        for ii=1:nn
            Cij(ii,:)=gammat(n+(1:nn)-ii+1);
        end
        Cia=gammat((nn:-1:1)+n+1);
        l = transpose(Cij\Cia);
        S = g(1)-l*Cia;
        L{u}(i,:)=simplify([zeros(1,i-nn-1) -l 1 zeros(1,n-i)]/sqrt(S));
    end
end



r=2.5;
Ct =cell(size(f,1),numel(L));
true=cell(size(f,1),1);
for j=1:size(f,1)
    true{j}=g(1:end-1);
    for i=2:numel(g)
        if j==3
            hh = min((i-1)*1,r/f{j,2});
        else
            hh = (i-1);
        end
        true{j} = subs(true{j},g(i),f{j,1}(hh, r/f{j,2}));
    end

    for u=1:numel(neigh)
        Lt = L{u};
        for i=2:numel(g)
            if j==3
            hh = min((i-1)*1,r/f{j,2});
        else
            hh = (i-1);
        end
            Lt = subs(Lt,g(i),f{j,1}(hh, r/f{j,2}));
        end
        Lt = eval(Lt);
        Ct{j,u} = Lt \ transpose(inv(Lt));
    end
end

figure(1);
for j=1:size(f,1)
    subplot(3,2,j); hold off 
    for u=1:numel(neigh)
        plot(0:n-1,Ct{j,u}(:,1)); hold on;
    end
    plot(0:n-1,true{j},'--k')
    axis tight; ylim([0 1]); xlabel(f{j,3}); ylabel('Covariogram')
    temp= textscan(num2str(neigh),'%s');
    legend(temp{1},'Theorical Equation')
    
    % set(gca,'xtick',[]); set(gca,'ytick',[])
end

%% Define path

n=9;
neigh = 2;
g = [1; sym('g%d', [n 1])];
syms h rr
f = {
    symfun(exp(-h/rr)                                      ,[h,rr]), 1, 'exponential';
    symfun(exp(-(h/rr).^2)                                 ,[h,rr]), 2/(pi^.5), 'gaussian';
    symfun(1 - 1.5*h/rr + .5*(h/rr).^3                     ,[h,rr]), 8/3, 'spherical';
    symfun(exp(-sqrt(h/rr))                                ,[h,rr]), 1/2, 'stable';
    symfun(1./(1+h/rr)                                     ,[h,rr]), 1/4, 'hyperbolic';
    symfun(sin(h/rr)./(h/rr)                               ,[h,rr]) 2/3, 'cardinal sine';
    };

% Mid-point
path=[1 9 5 3 7 2 4 6 8];
% path=[1 3 2];
sim=false(n,1);
x=sort(path);
for u=1:numel(neigh)
    L{u} = sym('L%d%d', [n n]);
    for i=1:n
        x_sim = path(i);
        x_neigh = x(sim);
        [~,idx] = sort(abs(x_neigh-x_sim));
        x_neigh = x_neigh(idx);
        x_neigh = x_neigh(1:min(sum(sim),neigh(u)));
        Cij = g(pdist2(x_neigh',x_neigh')+1);
        Cia = g(pdist2(x_neigh',x_sim')+1);
        l = transpose(Cij\Cia);
        S = g(1)-l*Cia;
        L{u}(path(i),:)=0;
        L{u}(path(i),path(i))=1/sqrt(S);
        L{u}(path(i),x_neigh)=-l/sqrt(S);
        sim(path(i))=true;
    end
end

C = L{1} \ transpose(inv(L{1}));

%% All possible path
n=9;
neigh = 2;
rr = 2;
vario = {'exponential', 'gaussian','spherical', 'hyperbolic', 'cardinal sine'};
parm.gen.covar(1).range = 2;
parm.gen.covar(1).azimuth = [];
parm.gen.covar(1).c0 = 1;

x=(1:n)'; %distance position
for v=1:numel(vario)
    parm.gen.covar.model = vario{v};
    parm = kriginginitiaite(parm);
    Ct{v} = covardm_perso(x,x,parm.k.covar);
    g(v,:) = Ct{v}(1,:);
end

path = perms(x);

tic
for v=1:numel(vario)
    parfor j=1:size(path,1)
        for u=1:numel(neigh)
            L = nan(n, n);
            sim=false(n,1);
            for i=1:n
                x_sim = path(j,i);
                x_neigh = x(sim);
                [~,idx] = sort(abs(x_neigh-x_sim));
                x_neigh = x_neigh(idx);
                x_neigh = x_neigh(1:min(sum(sim),neigh(u)));
                if numel(x_neigh)==0
                    l=0;
                    S=g(v,1);
                else
                    Cij = reshape(g(v,pdist2(x_neigh,x_neigh)+1),numel(x_neigh),numel(x_neigh));
                    Cia = g(v,pdist2(x_neigh,x_sim)+1)';
                    l = transpose(Cij\Cia);
                    S = g(v,1)-l*Cia;
                end
                L(path(j,i),:)=0;
                L(path(j,i),path(j,i))=1/sqrt(S);
                L(path(j,i),x_neigh)=-l/sqrt(S);
                sim(path(j,i))=true;
            end
            C{j} = L \ transpose(inv(L));
        end
    end
    CC{v}=C;
    toc
end

parfor j=1:size(path,1)
    fne(v,j)=sum((C{v}{j}(:)-Ct{v}(:)).^2);
end

[fne_s, idx] = sort(fne);
path_s =path(idx,:);

histogram(fne)