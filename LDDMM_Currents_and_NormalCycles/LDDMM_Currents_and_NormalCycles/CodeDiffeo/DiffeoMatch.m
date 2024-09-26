function DF = DiffeoMatch(H,gammaR)

% This implements the LDDMM matching functional and its gradient :
% See script_match.m for an example of use
% The inputs are
% H : output to a call to HamiltSys function, which contains settings and equations for
% the model of deformations
% gammaR : a non negative scalar giving the weight of the energy term

InitPos = [];
    function q0 = GetInitPos
        q0 = InitPos;
    end
DF.GetInitPos = @GetInitPos;

Target = {};
TargetTime = [];
TargetRange = {};
ntargets = 0;
if nargin<2
    gammaR = 0;
end

    function [n,d] = GetDims
        % dimensions of the data set up in the deformation model
        % n is the number of points and d the dimension
        [n,d] = H.GetDims();
    end
DF.GetDims = @GetDims;

    function SetDims(n,d)
        % dimensions of the data set up in the deformation model
        % n is the number of points and d the dimension
        H.SetDims(n,d);
    end
DF.SetDims = @SetDims;

    function SortTargets
        % internal function to sort targets according to their time
        [TargetTime,ind] = sort(TargetTime);
        Target = Target(ind);
        TargetRange = TargetRange(ind);
        ntargets = length(Target);
    end

    function AddTarget(tgt,t,tgt_range)
        % add target tgt with time parameter t to the functional
        if nargin<2
            t = 1;
        end
        Target{end+1} = tgt;
        TargetTime(end+1) = t;
        % get initial source points and test if they match the initial 
        % source points of one of the previous targets
        tgt_npoints = size(tgt.InitPos,1);
        n = size(InitPos,1);
        
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
            H.SetDims(n+tgt_npoints,size(InitPos,2));
        end
        SortTargets();
    end
DF.AddTarget = @AddTarget;

    function [target,targettime,targetrange] = GetTargets()
        target = Target;
        targettime = TargetTime;
        targetrange = TargetRange;
    end
DF.GetTargets = @GetTargets;

    function J = Funct(p,q)
        % computes the functional J as a function of the inital momentum  p
        % and initial positions q. J(p,q) is computed by summing the energy 
        % term and the target functionals
        t = 0;
        J = 0;
        % first computes energy term
        if gammaR
            J = H.Energy(p,q,gammaR);
        end
        % loop through targets which are sorted according to time
        for k = 1:ntargets
            if TargetTime(k)>t
                % perform shooting to move to the next target time
                [p,q] = H.Shoot(p,q,t,TargetTime(k));
                t = TargetTime(k);
            end
            J = J + Target{k}.Eval(q(TargetRange{k},:));
        end
    end
DF.Funct = @Funct;

    function [u,v] = GradFunct(p,q)
        % gradient of functional J as a function of initial momentum p and
        % positions q.
        t = 1;
        % perform shooting to move from time 0 to time 1
        [p,q] = H.Shoot(p,q);
        % gradient of the energy term
        if gammaR
            [u,v] = H.GradEnergy(p,q,gammaR);
        else
            u = zeros(size(p));
            v = u;
        end
        % loop through targets in reverse time ordering
        for k = ntargets:-1:1
            if TargetTime(k)<t
                % shoot backward the gradient to the previous target time
                [~,q,u,v] = H.ShootGrad_Precomp(u,v,t,TargetTime(k));
                t = TargetTime(k);
            end
            % add the current target gradient
            v(TargetRange{k},:) = v(TargetRange{k},:) + Target{k}.Grad(q(TargetRange{k},:));
        end
        if t>0
            % finally shoot backward the gradient to time 0
            [~,~,u,v] = H.ShootGrad_Precomp(u,v,t,0);
        end
    end
DF.GradFunct = @GradFunct;

    function [fSetTime,fToggleGrid] = Plot(p0,q0,ctag,withgrid)
        % function for plotting the points trajectories and targets
        if nargin<4
            withgrid = 1;
        end
        if nargin<3
            ctag = 'b';
        end
        hold on
        [fSetTimeH,funGrid] = H.PlotShoot(p0,q0,ctag,withgrid);
        axis off
        marker = {'s','x','+','*','.','d','^','v','<','>','p','h'};
        marker(end:ntargets) = {'x'};
        t = 0; p = p0; q = q0;
        for k = 1:ntargets
            if TargetTime(k)>t
                % perform shooting to move to the next target time
                [p,q] = H.Shoot(p,q,t,TargetTime(k));
                t = TargetTime(k);
            end
            h{k} = Target{k}.Plot(q(TargetRange{k},:),[marker{k},ctag]);
        end
        shg
        function ToggleGrid
            v = axis;
            funGrid();
            axis(v);
        end
        fToggleGrid = @ToggleGrid;
        function SetTime(t)
            fSetTimeH(t);
            tc = 0; p = p0; q = q0;
            for k = 1:ntargets
                if TargetTime(k)>tc
                    % perform shooting to move to the next target time
                    tn = min(t,TargetTime(k));
                    if tn>tc
                        [p,q] = H.Shoot(p,q,tc,tn);
                    end
                    tc = tn;
                end
                delete(h{k})
                h{k} = Target{k}.Plot(q(TargetRange{k},:),[marker{k},ctag]);
            end
        end
        fSetTime = @SetTime;
    end
DF.Plot = @Plot;

end