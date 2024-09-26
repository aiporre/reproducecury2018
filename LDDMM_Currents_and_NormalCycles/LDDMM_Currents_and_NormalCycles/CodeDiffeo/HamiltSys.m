function H = HamiltSys(ker,n,d,solver,odesettings)

% Implements geodesic equations and gradient geodesic equations
% for a given kernel
% n is number of points
% d is dimension

% default dimensions (to be updated afterwards via SetDims function)
if nargin<3
    d = 1;
end
if nargin<2
    n = 1;
end

    function SetKer(ker_)
        ker = ker_;
    end
H.SetKer = @SetKer;

    function SetDims(n_,d_)
        n = n_;
        d = d_;
        SetRanges();
    end
H.SetDims = @SetDims;

    function [n_,d_] = GetDims
        n_ = n;
        d_ = d;
    end
H.GetDims = @GetDims;

% solver settings
if nargin<4
    solver = @ode45;
end
if nargin<5
    odesettings = odeset();%odeset('MaxStep',1);%odeset('RelTol',1e-6,'AbsTol',1e-9);%
end


%%%%%%%%%%%%%%%%%%%%%%%%
%  Geodesic equations  %
%%%%%%%%%%%%%%%%%%%%%%%%

    function [pt,qt,T,P,Q,sol] = Shoot(ps,qs,s,t)
        % integrate geodesic equations from time s to time t
        if nargin<3
            s = 0;
            t = 1;
        end
        SolShoot = solver(@GeodEq,[s,t],join(ps,qs),odesettings);
        T = SolShoot.x;
        
        % debugging : check that hamiltonian is constant over time
        %clf;plot((sum(SolShoot.y(1:end/2,1:end-1).*diff(SolShoot.y(end/2+1:end,:),1,2))./diff(T)));pause

        [P,Q] = split(SolShoot.y);
        pt = P(:,:,end);
        qt = Q(:,:,end);
        if nargout==6
            sol = SolShoot;
        end
    end
H.Shoot = @Shoot;

    function pqdot = GeodEq(~,pq)
        [p,q] = split(pq);
        ker.Precomp(q);
        pdot = -.5 * ker.Grad_Precomp(p,p);
        qdot = ker.Eval_Precomp(p);
        pqdot = join(pdot,qdot);
    end




%%%%%%%%%%%%%%%%%%%%
%  Flow equations  %
%%%%%%%%%%%%%%%%%%%%

    function [xt,X] = Flow(x,ps,qs,s,t,nt)
        % integrate flow equations from time s to time t
        % ps, qs : initial positions and momenta (n*d matrices)
        % x : points to be flowed (nx*d matrix)
        % nt : number of time points for the full trajectories in X
        nx = size(x,1);
        xrange = 2*n*d+(1:nx*d);
        if nargin<6
            nt = 2;
        end
        if nargin<5
            s = 0;
            t = 1;
        end
        tspan = linspace(s,t,nt);
        [~,PQX] = solver(@FlowEq,tspan,join_flow(ps,qs,x),odesettings);
        [~,~,X] = split_flow(PQX');
        xt = X(:,:,end);
    end
H.Flow = @Flow;

    function pqxdot = FlowEq(~,pqx)
        [p,q,x] = split_flow(pqx);
        ker.Precomp(q);
        pdot = -.5 * ker.Grad_Precomp(p,p);
        qdot = ker.Eval_Precomp(p);
        xdot = ker.Eval(x,q,p);
        pqxdot = join_flow(pdot,qdot,xdot);
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Gradient Geodesic equations  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function [pt,qt,ut,vt] = ShootGrad(p,q,u,v,s,t)
        [~,tPQUV] = solver(@GradEq,[s,t],join2(p,q,u,v),odesettings);
        [pt,qt,ut,vt] = split2(tPQUV(end,:)');
    end
H.ShootGrad = @ShootGrad;

    function pquvdot = GradEq(~,pquv)
        [p,q,u,v] = split2(pquv);
        ker.Precomp(q);
        pdot = -.5 * ker.Grad_Precomp(p,p);
        qdot = ker.Eval_Precomp(p);
        udot = ker.Diff_Precomp(p,u) - ker.Eval_Precomp(v);
        vdot = .5*ker.Hess_Precomp(p,u) - ker.Grad_Precomp(p,v);
        pquvdot = join2(pdot,qdot,udot,vdot);
    end


% Precomputed version (to avoid recomputing p and q when the geodesic
% equations have already been solved

SolShoot = []; % This variable is used to store the result of Shoot

    function [pt,qt,ut,vt] = ShootGrad_Precomp(u,v,s,t)
        [~,tUV] = solver(@GradEq_Precomp,[s,t],join(u,v),odesettings);
        [ut,vt] = split(tUV(end,:)');
        [pt,qt] = split(deval(SolShoot,t));
    end
H.ShootGrad_Precomp = @ShootGrad_Precomp;

    function uvdot = GradEq_Precomp(t,uv)
        [u,v] = split(uv);
        [p,q] = split(deval(SolShoot,t));
        ker.Precomp(q);
        udot = ker.Diff_Precomp(p,u) - ker.Eval_Precomp(v);
        vdot = .5*ker.Hess_Precomp(p,u) - ker.Grad_Precomp(p,v);
        uvdot = join(udot,vdot);
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Energy and gradient of energy  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function J = Energy(p,q,weight)
        J = weight * sum(sum(ker.Eval(q,q,p).*p));
    end
H.Energy = @Energy;

    function [u,v] = GradEnergy(p,q,weight)
        u = 2 * weight * ker.Eval(q,q,p);
        v = weight * ker.Grad(q,p,p);
    end
H.GradEnergy = @GradEnergy;



%%%%%%%%%%%%%%
%  Plotting  %
%%%%%%%%%%%%%%

    function [funSetTime,fToggleGrid] = PlotShoot(p0,q0,ctag,withgrid,withtrajs)
        if nargin<5
            withtrajs = 1;
        end
        if nargin<4
            withgrid = 1;
        end
        if nargin<3
            ctag = 'b';
        end
        hold on
        [~,~,~,~,~,SolShoot] = Shoot(p0,q0);
        nT = 100;
        PQ = deval(SolShoot,linspace(0,1,nT));
        Q = reshape(PQ(1+d*n:end,:),n,d,nT);
        if withgrid
            [funSetTimeGrid,funGrid] = PlotGrid(p0,q0,ctag);
        end
        if withtrajs
            h = PlotTrajs(Q,ctag);
        else
            h = [];
        end
        function ToggleGrid
            if ~exist('funGrid','var')
                [funSetTimeGrid,funGrid] = PlotGrid(p0,q0,ctag);
                withgrid = 1;
            else
                funGrid();
                withgrid = ~withgrid;
            end
        end
        fToggleGrid = @ToggleGrid;
        function SetTime(t)
            if withgrid
                funSetTimeGrid(t);
            end
            delete(h)
            h = PlotTrajs(Q,ctag,0,t);
        end
        funSetTime = @SetTime;
        axis equal
    end
H.PlotShoot = @PlotShoot;

    function [funSetTime,fToggleGrid] = PlotGrid(p0,q0,ctag)
        withgrid = 1;
        if nargin<3
            ctag = 'g';
        end
        switch d
            case 2
                ng = 50;
            case 3
                ng = 10;
        end
        % set dimensions of enclosing box
        ming = min(q0) - max(max(q0)-min(q0))*0.5;
        maxg = max(q0) + max(max(q0)-min(q0))*0.5;
        sz = maxg-ming;
        dg = max(sz)/ng;
        h = [];
        hold on
        nT = 100;
        switch d
            case 2
                [X,Y] = meshgrid(ming(1):dg:maxg(1),ming(2):dg:maxg(2));
                sX = size(X);
                [phiXY,PhiXY] = Flow([X(:),Y(:)],p0,q0,0,1,nT);
                phiX = reshape(phiXY(:,1),sX);
                phiY = reshape(phiXY(:,2),sX);
                h = [h;plot(phiX,phiY,ctag)];
                h = [h;plot(phiX',phiY',ctag)];
                rgb = get(h(1),'Color');
                set(h,'Color',rgb+(1-rgb)*.5);
            case 3
                [X,Y,Z] = meshgrid(ming(1):dg:maxg(1),ming(2):dg:maxg(2),ming(3):dg:maxg(3));
                sX = size(X);
                [phiXYZ,PhiXYZ] = Flow([X(:),Y(:),Z(:)],p0,q0,0,1,nT);
                phiX = reshape(phiXYZ(:,1),sX);
                phiY = reshape(phiXYZ(:,2),sX);
                phiZ = reshape(phiXYZ(:,3),sX);
                for k = 1:3
                    h = [h;plot3(reshape(phiX,sX(1),sX(2)*sX(3)),...
                        reshape(phiY,sX(1),sX(2)*sX(3)),...
                        reshape(phiZ,sX(1),sX(2)*sX(3)),ctag)];
                    phiX = shiftdim(phiX,1);
                    phiY = shiftdim(phiY,1);
                    phiZ = shiftdim(phiZ,1);
                    sX = sX([2,3,1]);
                end
                rgb = get(h(1),'Color');
                set(h,'Color',rgb+(1-rgb)*.75);
        end
        function SetTime(t)
            switch d
                case 2
                    phiXY = PhiXY(:,:,round(t*(nT-1)+1));
                    phiX = reshape(phiXY(:,1),sX);
                    phiY = reshape(phiXY(:,2),sX);
                    delete(h)
                    h = [];
                    h = [h;plot(phiX,phiY,ctag)];
                    h = [h;plot(phiX',phiY',ctag)];
                    rgb = get(h(1),'Color');
                    set(h,'Color',rgb+(1-rgb)*.5);
                case 3
                    phiXYZ = PhiXYZ(:,:,round(t*(nT-1)+1));
                    phiX = reshape(phiXYZ(:,1),sX);
                    phiY = reshape(phiXYZ(:,2),sX);
                    phiZ = reshape(phiXYZ(:,3),sX);
                    delete(h)
                    h = [];
                    for k = 1:3
                        h = [h;plot3(reshape(phiX,sX(1),sX(2)*sX(3)),...
                            reshape(phiY,sX(1),sX(2)*sX(3)),...
                            reshape(phiZ,sX(1),sX(2)*sX(3)),ctag)];
                        %set(h,'Color',rgbcolor);
                        phiX = shiftdim(phiX,1);
                        phiY = shiftdim(phiY,1);
                        phiZ = shiftdim(phiZ,1);
                        sX = sX([2,3,1]);
                    end
                    rgb = get(h(1),'Color');
                    set(h,'Color',rgb+(1-rgb)*.75);
            end
            uistack(h,'bottom');
        end
        funSetTime = @SetTime;
        function ToggleGrid
            withgrid = ~withgrid;
            if withgrid
                set(h,'Visible','on');
            else
                set(h,'Visible','off');
            end
        end
        fToggleGrid = @ToggleGrid;
    end
H.PlotGrid = @PlotGrid;

    function h = PlotTrajs(X,tag,marker,t)
        % plot trajectories
        if nargin<4
            t = 1;
        end
        if nargin<3
            marker  = 0;
        end
        if nargin<2
            tag = '';
        end
        hold on
        x0 = X(:,:,1);
        nT = size(X,3);
        rT = 1:round(t*(nT-1)+1);
        switch size(X,2)
            case 2
                h = plot(squeeze(X(:,1,rT))',squeeze(X(:,2,rT))',tag,'LineWidth',3);
                if marker
                    h = [h;plot(x0(:,1),x0(:,2),marker)];
                end
            case 3
                h = plot3(squeeze(X(:,1,rT))',squeeze(X(:,2,rT))',squeeze(X(:,3,rT))',tag,'LineWidth',3);
                if marker
                    h = [h;plot3(x0(:,1),x0(:,2),x0(:,3),marker)];
                end
        end
        axis equal
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  functions for reshaping data  %
%  for use with solvers          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prange = [];
qrange = [];
urange = [];
vrange = [];
nx = [];
xrange = [];
SetRanges();

    function SetRanges()
        prange = 1:n*d;
        qrange = n*d+prange;
        urange = n*d+qrange;
        vrange = n*d+urange;
    end

    function [P,Q] = split(PQ)
        % PQ : (2nd)*nT => P,Q : n*d*nT
        nT = size(PQ,2);
        P = reshape(PQ(prange,:),n,d,nT);
        Q = reshape(PQ(qrange,:),n,d,nT);
    end

    function PQX = join_flow(P,Q,X)
        % P,Q : n*d*nT, X : nx*d*nT => PQX : (2nd+nxd)*nT
        nT = size(P,3);
        PQX = [reshape(P,n*d,nT);reshape(Q,n*d,nT);reshape(X,nx*d,nT)];
    end

    function [P,Q,X] = split_flow(PQX)
        % PQX : (2nd+md)*nT => P,Q : n*d*nT, X : m*d*nT
        nT = size(PQX,2);
        P = reshape(PQX(prange,:),n,d,nT);
        Q = reshape(PQX(qrange,:),n,d,nT);
        X = reshape(PQX(xrange,:),nx,d,nT);
    end

    function PQ = join(P,Q)
        % P,Q : n*d*nT => PQ : (2nd)*nT
        nT = size(P,3);
        PQ = [reshape(P,n*d,nT);reshape(Q,n*d,nT)];
    end

    function [P,Q,U,V] = split2(PQUV)
        % PQUV : (4nd)*nT => P,Q,U,V : n*d*nT
        nT = size(PQUV,2);
        P = reshape(PQUV(prange,:),n,d,nT);
        Q = reshape(PQUV(qrange,:),n,d,nT);
        U = reshape(PQUV(urange,:),n,d,nT);
        V = reshape(PQUV(vrange,:),n,d,nT);
    end

    function PQUV = join2(P,Q,U,V)
        % P,Q,U,V : n*d*nT => PQUV : (4nd)*nT
        nT = size(P,3);
        PQUV = [reshape(P,n*d,nT);reshape(Q,n*d,nT);reshape(U,n*d,nT);reshape(V,n*d,nT)];
    end


end