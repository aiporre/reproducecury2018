function tgt = TargetMeasures(x,y,ker,wx,wy)

% Measure Matching target : kernel squared distance between a set of
% moving points q and the set of target points
% input x is a nx*d matrix containing the coordinates of the source points
% input y is a ny*d matrix containing the coordinates of the target points
% input ker is the kernel used for computing the distance
% input wx is a vector of length nx of nonnegative numbers giving weights for
% each point in x
% input wy is a vector of length y of nonnegative numbers giving weights for
% each point in y

tgt.InitPos = x;

nx = size(x,1);
ny = size(y,1);
nz = nx + ny;

if nargin < 4
    wx = ones(nx,1)/nx;
    wy = ones(ny,1)/ny;
end

wz = [wx;-wy];

% computes the measure matching functional
    function A = Eval(q)
        z = [q;y];
        A = wz'*ker.Eval(z,z,wz);
    end
tgt.Eval = @Eval;

% computes the gradient of the measure matching functional
    function G = Grad(q)
        z = [q;y];
        G = zeros(size(q));
        Gz = ker.Grad(z,wz,wz);
        G = Gz(1:nx,:);
    end
tgt.Grad = @Grad;

% plotting function
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

