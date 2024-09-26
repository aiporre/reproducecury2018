function DF = SmallDefMatch(ker,gammaR)

% This implements small deformation matching functional and its gradient :
% See script_match_smalldef.m for an example of use
% The inputs are
% ker : kernel for the model of deformations
% gammaR : a non negative scalar giving the weight of the energy term

InitPos = [];
    function q0 = GetInitPos
        q0 = InitPos;
    end
DF.GetInitPos = @GetInitPos;

Target = {};
TargetRange = {};
ntargets = 0;
if nargin<2
    gammaR = 0;
end

n = 0;
d = 0;

    function [n_,d_] = GetDims
        % dimensions of the data set up in the deformation model
        % n is the number of points and d the dimension
        n_ = n;
        d_ = d;
    end
DF.GetDims = @GetDims;

    function SetDims(n_,d_)
        % dimensions of the data set up in the deformation model
        % n is the number of points and d the dimension
        n = n_;
        d = d_;
    end
DF.SetDims = @SetDims;

    function AddTarget(tgt,tgt_range)
        % add target tgt to the functional
        Target{end+1} = tgt;
        % get initial source points and test if they match the initial 
        % source points of one of the previous targets
        tgt_npoints = size(tgt.InitPos,1);
        test = 1;
        if nargin==3
            TargetRange{end+1} = tgt_range;
            test = 0;
        else
            for k = 1:ntargets
                if tgt_npoints==length(TargetRange{k}) && all(all(InitPos(TargetRange{k},:)==tgt.InitPos))
                    TargetRange{end+1} = TargetRange{k};
                    test = 0;
                    break;
                end
            end
        end
        if test
            TargetRange{end+1} = n+1:n+tgt_npoints;
            InitPos = [InitPos;tgt.InitPos];
        end
        [n,d] = size(InitPos);
        ntargets = ntargets + 1;
    end
DF.AddTarget = @AddTarget;

    function [target,targetrange] = GetTargets()
        target = Target;
        targetrange = TargetRange;
    end
DF.GetTargets = @GetTargets;

    function J = Funct(p,q)
        % computes the functional J as a function of the inital momentum  p
        % and initial positions q. J(p,q) is computed by summing the energy 
        % term and the target functionals
        J = 0;
        qdot = ker.Eval(q,q,p);
        % first computes energy term
        if gammaR
            J = gammaR * sum(qdot(:).*p(:));
        end
        q = q + qdot;
        % loop through targets
        for k = 1:ntargets
            J = J + Target{k}.Eval(q(TargetRange{k},:));
        end
    end
DF.Funct = @Funct;

    function [Gp,Gq] = GradFunct(p,q)
        % gradient of functional J as a function of initial momentum p and
        % positions q.
        ker.Precomp(q);
        q1 = q + ker.Eval_Precomp(p);
        G = zeros(size(p));
        for k = 1:ntargets
            % add the current target gradient
            G(TargetRange{k},:) = G(TargetRange{k},:) + Target{k}.Grad(q1(TargetRange{k},:));
        end
        Gp = ker.Eval_Precomp(2*gammaR*p+G);
        if nargout>1
            tmp = gammaR*p+G;
            Gq = G + ker.Grad_Precomp(tmp,p);
        end
    end
DF.GradFunct = @GradFunct;

    function [fSetTime,fToggleGrid] = Plot(p,q,ctag,withgrid)
        % function for plotting the points trajectories and targets
        if nargin<4
            withgrid = 1;
        end
        t = 1;
        h = [];
        if nargin<3
            ctag = 'b';
        end
        PlotGrid = [];
        funToggleGrid = [];
        if withgrid
            [PlotGrid,funToggleGrid] = SetPlotGrid(p,q,ctag);
        end
        function SetTime(t_)
            t = t_;
            DrawNow;
        end
        fSetTime = @SetTime;
        q1 = q + ker.Eval(q,q,p);
        function DrawNow()
            delete(h)
            h = [];
            qt = (1-t)*q + t*q1;
            hold on
            if withgrid
                PlotGrid(t);
            end
            switch d
                case 2
                    h = [h;plot([q(:,1),qt(:,1)]',[q(:,2),qt(:,2)]',ctag,'LineWidth',3)];
                case 3
                    h = [h;plot3([q(:,1),qt(:,1)]',[q(:,2),qt(:,2)]',[q(:,3),qt(:,3)]',ctag,'LineWidth',3)];
            end
            axis off
            marker = {'s','x','+','*','.','d','^','v','<','>','p','h'};
            marker(end:ntargets) = {'x'};
            for k = 1:ntargets
                h = [h;Target{k}.Plot(qt(TargetRange{k},:),[marker{k},ctag])];
            end
        end
        DrawNow;
        function ToggleGrid
            if isempty('PlotGrid')
                [PlotGrid,funToggleGrid] = SetPlotGrid(p,q,ctag);
                PlotGrid();
                withgrid = 1;
            else
                funToggleGrid();
                withgrid = ~withgrid;
            end
        end
        fToggleGrid = @ToggleGrid;
    end
DF.Plot = @Plot;

    function [fDrawNow,fToggleGrid] = SetPlotGrid(p0,q0,ctag)
        withgrid = 0;
        h = [];
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
        switch d
            case 2
                [X,Y] = meshgrid(ming(1):dg:maxg(1),ming(2):dg:maxg(2));
                sX = size(X);
                XY = [X(:),Y(:)];
                phiXY = XY + ker.Eval(XY,q0,p0);
            case 3
                [X,Y,Z] = meshgrid(ming(1):dg:maxg(1),ming(2):dg:maxg(2),ming(3):dg:maxg(3));
                sX = size(X);
                XYZ = [X(:),Y(:),Z(:)];
                phiXYZ = XYZ + ker.Eval(XYZ,q0,p0);
        end
        function DrawNow(t)
            withgrid = 1;
            delete(h)
            h = [];
            hold on
            switch d
                case 2
                    phiXYt = (1-t)*XY + t*phiXY;
                    phiX = reshape(phiXYt(:,1),sX);
                    phiY = reshape(phiXYt(:,2),sX);
                    h = [h;plot(phiX,phiY,ctag)];
                    h = [h;plot(phiX',phiY',ctag)];
                    rgb = get(h(1),'Color');
                    set(h,'Color',rgb+(1-rgb)*.5);
                case 3
                    phiXYZt = (1-t)*XYZ + t*phiXYZ;
                    phiX = reshape(phiXYZt(:,1),sX);
                    phiY = reshape(phiXYZt(:,2),sX);
                    phiZ = reshape(phiXYZt(:,3),sX);
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
        end
        fDrawNow = @DrawNow;
        function ToggleGrid
            withgrid = ~withgrid;
            if withgrid
                v = axis;
                set(h,'Visible','on');
                axis(v)
            else
                v = axis;
                set(h,'Visible','off');
                axis(v);
            end
        end
        fToggleGrid = @ToggleGrid;
    end
DF.PlotGrid = @PlotGrid;


end