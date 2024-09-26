function OD = OptimDiffeo(DM,maxiter)

% This contains the optimization functions to perform matching or template
% estimation. It is a kind of wrapper for the fminlbfgs function
% input is DM the output to a call to DiffeoMatch
% or a cell array containg several calls to DiffeoMatch
% see script_match.m and script_template.m for examples of use.

if nargin < 3
    debug_testgrad = 0;
end

if nargin < 2
    maxiter = inf;
end

if ~iscell(DM)
    % if input is not a cell array,we convert to cell array of length 1
    tmp{1} = DM;
    DM = tmp;
    clear tmp
end

% n is number of points, d is dimension
[n,d] = DM{1}.GetDims();

% N is number of Deformations we need to optimize
% N=1 means we optimze only one deformation (matching)
% N>1 means we will do template estimation
N = length(DM);

q0 = [];

% options structure for the fminlbfgs function
lbfgs_options = struct('GoalsExactAchieve',1 ...
    ,'GradConstr',0 ...
    ,'HessUpdate','lbfgs' ...%'steepdesc'...%
    ,'StoreN',20 ...
    ,'GradObj','on' ...
    ,'MaxFunEvals',inf ...
    ,'MaxIter',maxiter ... 
    ,'Display','iter' ... %'bar' ... % 'off' ... %  'notify'
    ,'TolX',1e-5 ...
    ,'TolFun',1e-5 ...
    );
%    ,'OutputFcn',@Plot_fminlbfgs ...
%                  ,'rho',0 ...
%                  ,'sigma',inf ...

    function SetOptimOption(fieldname,fieldval)
        if isnumeric(fieldval)
            fieldval = num2str(fieldval);
        else
            fieldval = ['''',fieldval,''''];
        end
        eval(['lbfgs_options.',fieldname,' = ',fieldval,';']);
    end

OD.SetOptimOption = @SetOptimOption;

% functions for converting several momentum vectors to one single vector
% and conversely
    function p_v = cell2vect(p)
        range = 1:n*d;
        p_v = zeros(N*n*d,1);
        for k = 1:N
            p_v(range) = p{k}(:);
            range = range + n*d;
        end
    end

    function p = vect2cell(p_v)
        range = 1:n*d;
        p = cell(1,N);
        for k = 1:N
            p{k} = reshape(p_v(range),n,d);
            range = range + n*d;
        end
    end            

% performs the optimization by calling the fminlbfgs function. This
% optimizes both initial momentum vectors and the template points for
% template estimation
    function [p0,q0,J,exitflag,output,grad] = Optimize(p0,q0)
        [n,d] = DM{1}.GetDims();
        if N==1 && ~iscell(p0)
            tmp{1} = p0;
            p0 = tmp;
        end
        [pq0,J,exitflag,output,grad] = fminlbfgs(@fun_fminlbfgs,[cell2vect(p0);q0(:)],lbfgs_options);
        p0 = vect2cell(pq0(1:N*d*n));
        if N==1
            p0 = p0{1};
        end
        q0 = reshape(pq0(1+d*n*N:end),n,d);
    end
OD.Optimize = @Optimize;

% performs the optimization by calling the fminlbfgs function. This
% optimizes only the initial momentum vectors for matching
    function [p0,J,exitflag,output,grad] = Optimize_onlymom(p0,q0_)
        [n,d] = DM{1}.GetDims();
        if N==1 && ~iscell(p0)
            tmp{1} = p0;
            p0 = tmp;
        end
        q0 = q0_;
        [n,d] = size(q0);
        for k=1:length(DM)
            DM{k}.SetDims(n,d);
        end
        [p0,J,exitflag,output,grad] = fminlbfgs(@fun_onlymom_fminlbfgs,cell2vect(p0),lbfgs_options);
        %[p0,J,exitflag,output,grad] = fmingradesc(@fun_onlymom_fminlbfgs,cell2vect(p0),.001,10);
        p0 = vect2cell(p0);
        if N==1
            p0 = p0{1};
        end
    end
OD.Optimize_onlymom = @Optimize_onlymom;

% function to compute functional and gradient to be passed to the fminlbfgs
% routine. It extracts p0 and q0 from the single vector pq0, and then
% computes and collects all LDDMM functionals.
    function [J,G] = fun_fminlbfgs(pq0)
        p0 = vect2cell(pq0(1:N*n*d));
        q0 = reshape(pq0(1+d*n*N:end),n,d);
        J = 0;
        for k = 1:N
            J = J + DM{k}.Funct(p0{k},q0);
        end
        if nargout > 1
            Gp = cell(1,N);
            Gq = 0;
            for k = 1:N
                [Gp{k},Gqk] = DM{k}.GradFunct(p0{k},q0);
                Gq = Gq + Gqk;
            end
            G = [cell2vect(Gp);Gq(:)];
        end
    end

% Same as previously but with optimization on momentum only.
    function [J,G] = fun_onlymom_fminlbfgs(p0)
        % function to compute functional and gradient
        % to be used with fminlbfgs
        p0 = vect2cell(p0);
        J = 0;
        for k = 1:N
            J = J + DM{k}.Funct(p0{k},q0);
        end
        if nargout > 1
            Gp = cell(1,N);
            for k = 1:N
                Gp{k} = DM{k}.GradFunct(p0{k},q0);
            end
            G = cell2vect(Gp);
%             if 1
%                 %test gradient (for debugging)
%                 p0{1}(1:10,:)'
%                 p0=cell2vect(p0);
%                 lambda = [-100:-1,1:100]/10000;
%                 %lambda = [-1,1]/100000;
%                 testgrad = zeros(1,length(lambda));
%                 nG2 = norm(G)^2
% %                 for k=1:length(lambda)
% %                     testgrad(k)=(fun_onlymom_fminlbfgs(p0-lambda(k)*G)-J)/(-lambda(k)*nG2);
% %                 end
%                 for k=1:length(lambda)
%                     testgrad(k)=fun_onlymom_fminlbfgs(p0-lambda(k)*G);
%                 end
%                 clf
%                 plot(lambda,testgrad);
%                 hold on
%                 plot(0,J,'o')
%                 pause
%             end
        end
    end

% functions for plotting the results

    function ret = Plot_fminlbfgs(pq0,~,tag)
        % function to plot during iterations
        % to be used with fminlbfgs
        p0 = vect2cell(pq0(1:N*n*d));
        q0 = reshape(pq0(1+d*n*N:end),n,d);
        clf
        hold on
        for k = 1:N
            DM{k}.Plot(p0{k},q0);
        end
        title(tag)
        ret = 0;
    end

    function [funAnim,fToggleGrids] = Plot(p0,q0,withgrid)
        if nargin<3
            withgrid = 1;
        end
        if N==1 && ~iscell(p0)
            tmp{1} = p0;
            p0 = tmp;
        end
        ctag = {'b','r','g','k','c','y'};
        ctag(end:N) = {'m'};
        for k = 1:N
            [fSetTime{k},fToggleGrid{k}] = DM{k}.Plot(p0{k},q0,ctag{k},withgrid);
        end
        function ToggleGrids
            for k = 1:N
                fToggleGrid{k}();
            end
        end
        fToggleGrids = @ToggleGrids;
        function Anim
            v = axis;
            for t=linspace(0,1,30)
                for k = 1:N
                    fSetTime{k}(t);
                end
                axis(v)
                shg
            end
        end
        funAnim = @Anim;                                      
    end
OD.Plot = @Plot;

end
