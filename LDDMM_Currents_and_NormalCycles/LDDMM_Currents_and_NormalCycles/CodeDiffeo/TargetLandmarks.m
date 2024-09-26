function tgt = TargetLandmarks(x,y)

% This implements the Landmark Matching target : sum of squared distances
% of the points positions to their corresponding target points.
% input y is a n*d matrix containing the coordinates of the target points

tgt.InitPos = x;

% Evaluation function
    function A = Eval(q)
        A = sum((q(:)-y(:)).^2);
    end
tgt.Eval = @Eval;

% Gradient function
    function G = Grad(q)
        G = 2*(q-y);
    end
tgt.Grad = @Grad;

% plotting
    function h = Plot(q,tag)
        if nargin<1
            tag = 'x';
        end
        switch size(y,2)
            case 2
                h = [];%plot(q(:,1),q(:,2),tag);
                h = [h;plot(y(:,1),y(:,2),tag)];
            case 3
                h = [];%plot3(q(:,1),q(:,2),q(:,3),tag);
                h = [h;plot3(y(:,1),y(:,2),y(:,3),tag)];
        end
    end
tgt.Plot = @Plot;

end

