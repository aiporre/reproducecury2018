function tgt = TargetCurveCurr(S,T,ker)

% Curve Matching target : implements the currents metric between a moving
% curve and its target
% inputs S and T are structures corresponding to the source and target 
% curves. They must contain fields S.Vertices (resp. 
% T.Vertices) of size nvS*3 (resp. nvT*3) and S.Faces of size nfS*2 (resp. 
% T.Faces of size nfT*2)
% input ker is the kernel used for computing the distance

tgt.InitPos = S.Vertices;

% computes the vector dirac currents located at the centers of the 
% triangles of the target mesh
muT = Curve2Curr(T);

%sqn_muT = ScalProdCurr(muS,muS,ker);

% computes the currents functional
    function A = Eval(q)
        % updates the vertices coordinates of the moving mesh
        S.Vertices = q;
        % computes the vector dirac currents of the moving mesh
        muS = Curve2Curr(S);
        
        %A = ScalProdCurr(muS,muS,ker) - 2 * ScalProdCurr(muS,muT,ker) + sqn_muT;
        
        % computes the difference between the two sets of dirac currents
        mu = DiffCurr(muS,muT);
        % computes the squared kernel norm of the difference
        A = ScalProdCurr(mu,mu,ker);
    end
tgt.Eval = @Eval;

% computes the currents functional gradient
    function G = Grad(q)
        % updates the vertices coordinates of the moving mesh
        S.Vertices = q;
        % computes the vector dirac currents of the moving mesh
        muS = Curve2Curr(S);
        % computes the difference between the two sets of dirac currents
        mu = DiffCurr(muS,muT);
        % computes the gradient of the functional with respect to the
        % locations of the dirac currents (centers of the triangles)
        Gc = ker.Grad(mu.Points,mu.Vectors,mu.Vectors);
        Gc = Gc(1:size(muS.Points,1),:);
        % computes the gradient of the functional with respect to the
        % vectors of the dirac currents (normal vectors of the triangles)
        Gn = 2*ker.Eval(mu.Points,mu.Points,mu.Vectors);
        Gn = Gn(1:size(muS.Points,1),:);
        % transfer these two gradients to get the gradients at the vertices
        % of the moving mesh 
        G = TransCurve2Curr(S,Gc,Gn);
    end
tgt.Grad = @Grad;        

% plotting function
    function h = Plot(q,~)
        S.Vertices = q;
        h{1} = plotcurve(S,'g');
        h{2} = plotcurve(T,'r');
    end
tgt.Plot = @Plot;

% computes the set of Dirac currents (normals of the triangles located at
% the centers of the triangles) of a triangulated mesh.
    function mu = Curve2Curr(S)
        
        V1 = S.Vertices(S.Faces(:,1),:);
        V2 = S.Vertices(S.Faces(:,2),:);
        
        % compute centers of segments of curve S
        mu.Points = (1/2)*(V1 + V2);
        
        % compute tangent vectors of segments of curve S
        % tangent vectors are not normalized 
        % norm of each tangent vector gives length of segment
        mu.Vectors = V2-V1;
        % if weights are given then multiply tangent vectors by these weights
        if isfield(S,'Weights')
            mu.Vectors = mu.Vectors.*repmat(S.Weights(:),1,3);
        end
    end

% computes the difference between two sets of Dirac currents : this is just
% the concatenation of the two sets with negative signs for the
% second set of vectors
    function mud = DiffCurr(mu,nu)
        mud.Points = [mu.Points;nu.Points];
        mud.Vectors = [mu.Vectors;-nu.Vectors];
    end

% computes kernel scalar product between sets of dirac currents mu1 and mu2
    function s = ScalProdCurr(mu1,mu2,ker)
        s = sum(sum(mu1.Vectors.*ker.Eval(mu1.Points,mu2.Points,mu2.Vectors)));
    end

% function to transfer the gradients with respect to the centers and normals of the
% triangles to their vertices. It computes the transpose of the
% differential of the function Surf2Curr which computes the centers and
% normals from the vertices
    function G = TransCurve2Curr(S,Gc,Gt)
        
        G = zeros(size(S.Vertices));
        
        F1 = S.Faces(:,1);
        F2 = S.Faces(:,2);
        
        V1 = S.Vertices(F1,:);
        V2 = S.Vertices(F2,:);
        
        coef = (1/2)*Gc;
        nF = size(F1,1);
        for f = 1:nF
            coeff = coef(f,:);
            G(F1(f),:) = coeff;
            G(F2(f),:) = coeff;
        end;
        
        if isfield(S,'Weights')
            Gt = Gt.*repmat(S.Weights(:),1,3);
        end
        
        coef = Gt;
        for f = 1:nF
            coeff = coef(f,:);
            G(F1(f),:) = G(F1(f),:) - coeff;
            G(F2(f),:) = G(F2(f),:) + coeff;
        end
    end

end