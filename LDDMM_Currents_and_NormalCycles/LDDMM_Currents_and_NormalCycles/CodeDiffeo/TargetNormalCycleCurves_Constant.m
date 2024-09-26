function tgt = TargetNormalCycleCurves_Constant(C,S,sigmaW,CoeffRenorm)
%Kernel metric on normal cycles with constant normal kernel

if isfield(C,'x')
C.Vertices = C.x;
C.Faces = C.G;
S.Vertices = S.x;
S.Faces = S.G;
end

tgt.InitPos = C.Vertices;

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Spatial Kernel      %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function pConv = GaussKernelNC(sigma)
        % gaussian kernel
        % k(x,y)=exp(-|x-y|^2/sigma^2)
        % sigma : width of kernel
        oosigma2 = 1/sigma^2;
        function K = expo(x,y)
            m = size(x,1);
            n = size(y,1);
            K = sum(x.^2,2)*ones(1,n) + ones(m,1)*sum(y.^2,2)' - 2 * scal_prods(x,y);
            K = exp(-oosigma2*K);
        end
        pConv = @expo;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First derivtive Spatial Kernel %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function pConv = d1gauss(sigma)
        
        oosigma2 = 1/sigma^2;
        
        function K = d1expo(x,y)
            [m,d] = size(x);
            [n,d] = size(y);
            ymx = repmat(permute(y,[2,3,1]),[1,m,1]) - repmat(x',[1,1,n]);
            K = repmat((2*oosigma2)*exp(-oosigma2*sum(ymx.^2)),[d,1,1]).*ymx;
        end
        
        pConv = @d1expo;
        
    end

kerp = GaussKernelNC(sigmaW);
d1kerp = d1gauss(sigmaW);


%%%%%%%%%%%%%%%%%%%%%%%%%
%      Scalar prod      %
%%%%%%%%%%%%%%%%%%%%%%%%%

    function Scal = scal(C1,C2)
        V1 = C1.Vertices;
        V2 = C2.Vertices;
        n1 = size(V1,1);
        n2 = size(V2,1);
        F1 = C1.Faces;
        F2 = C2.Faces;
        
        kerPQ = kerp(V1,V2);
        E1 = ComputeEdgesC3D(C1);
        E2 = ComputeEdgesC3D(C2);
        NE1 = E1./repmat(sqrt(sum(E1.*E1,2)),1,3);
        NE2 = E2./repmat(sqrt(sum(E2.*E2,2)),1,3);
        
        A = zeros(n1,3);
        B = zeros(n2,3);
        
        for d = 1:3
            A(:,d) = accumarray(F1(:,2),NE1(:,d),[n1,1])- accumarray(F1(:,1),NE1(:,d),[n1,1]);
            B(:,d) = accumarray(F2(:,2),NE2(:,d),[n2,1])- accumarray(F2(:,1),NE2(:,d),[n2,1]);
        end
        AB = scal_prods(A,B);
        Scal = sigmaW^2*CoeffRenorm*pi^2/4*sum(sum(kerPQ.*AB));
    end

scalS = scal(S,S);
    function A = Eval(q)
        C.Vertices = q;
        A = scal(C,C) + scalS -2*scal(C,S);
    end
tgt.Eval = @Eval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradient of the scalar prod  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function grad = Grad(q)
        C.Vertices = q;
        C1 = C;
        C2 = S;
        %Extraction des elements utiles
        V1 = C1.Vertices;
        V2 = C2.Vertices;
        F1 = C1.Faces;
        F2 = C2.Faces;
        E1 = ComputeEdgesC3D(C1);
        E2 = ComputeEdgesC3D(C2);
        N = size(V1,1);
        M = size(V2,1);
        n = size(E1,1);
        m = size(E2,1);
        kerPQ = kerp(V1,V2);
        kerPP = kerp(V1,V1);
        normE1 = sqrt(sum(E1.*E1,2));
        NE1 = E1./repmat(normE1,1,3);
        normE2 = sqrt(sum(E2.*E2,2));
        NE2 = E2./repmat(normE2,1,3);
        d1kerpPQ = permute(d1kerp(V1,V2),[2,3,1]);
        d1kerpPP = permute(d1kerp(V1,V1),[2,3,1]);
        
        A = zeros(N,3);
        B = zeros(M,3);
        
        for d = 1:3
            A(:,d) = accumarray(F1(:,2),NE1(:,d),[N,1]) - accumarray(F1(:,1),NE1(:,d),[N,1]);
            B(:,d) = accumarray(F2(:,2),NE2(:,d),[M,1]) - accumarray(F2(:,1),NE2(:,d),[M,1]);
        end
        AA = scal_prods(A,A);
        AB = scal_prods(A,B);
        grad1PQ = sum(d1kerpPQ.*repmat(AB,1,1,3),2);  
        grad1PP = sum(d1kerpPP.*repmat(AA,1,1,3),2);
        
        Btemp = permute(repmat(B,1,1,n),[3,1,2]);
        Atemp = permute(repmat(A,1,1,n),[3,1,2]);
        normE1tempB = repmat(normE1,1,M,3);
        normE1tempA = repmat(normE1,1,N,3);
        %ABtemp = repmat(AB,1,1,3);
        NE1_Btemp = repmat(scal_prods(NE1,B),1,1,3);
        NE1_Atemp = repmat(scal_prods(NE1,A),1,1,3);
        %kerPQtemp = repmat(kerPQ,1,1,3);
        NE1tempB = permute(repmat(NE1,1,1,M),[1,3,2]);
        NE1tempA = permute(repmat(NE1,1,1,N),[1,3,2]);
        
        grad_edges_PP = (1./normE1tempA).*(Atemp - NE1_Atemp.*NE1tempA);
        grad_edges_PQ = (1./normE1tempB).*(Btemp - NE1_Btemp.*NE1tempB);
        grad_vertices_PQ = zeros(N,M,3);
        grad_vertices_PP = zeros(N,N,3);
        
        for d = 1:3       
            
            KerPPF12 = kerPP(F1(:,2),:);
            KerPPF11 = kerPP(F1(:,1),:);
            KerPQF12 = kerPQ(F1(:,2),:);
            KerPQF11 = kerPQ(F1(:,1),:);
            
            grad_edges_kerPP2 = KerPPF12.*grad_edges_PP(:,:,d);
            grad_edges_kerPP1 = KerPPF11.*grad_edges_PP(:,:,d);
            grad_edges_kerPQ2 = KerPQF12.*grad_edges_PQ(:,:,d);
            grad_edges_kerPQ1 = KerPQF11.*grad_edges_PQ(:,:,d);
            
            grad_vertices_PQ(:,:,d) = accummatrix(F1(:,2),grad_edges_kerPQ2,[N,M]) + accummatrix(F1(:,1),grad_edges_kerPQ1,[N,M]);
            grad_vertices_PQ(:,:,d) = grad_vertices_PQ(:,:,d) - accummatrix(F1(:,2),grad_edges_kerPQ1,[N,M]) - accummatrix(F1(:,1),grad_edges_kerPQ2,[N,M]);
            grad_vertices_PP(:,:,d) = accummatrix(F1(:,2),grad_edges_kerPP2,[N,N]) + accummatrix(F1(:,1),grad_edges_kerPP1,[N,N]);
            grad_vertices_PP(:,:,d) = grad_vertices_PP(:,:,d) - accummatrix(F1(:,2),grad_edges_kerPP1,[N,N]) - accummatrix(F1(:,1),grad_edges_kerPP2,[N,N]);
            
        end
        grad2PQ = sum(grad_vertices_PQ,2);
        grad2PP = sum(grad_vertices_PP,2);
        gradPQ = grad1PQ + grad2PQ;
        gradPP = grad1PP + grad2PP;
        grad = sigmaW^2*CoeffRenorm*2*pi^2/4*(squeeze(gradPP)-squeeze(gradPQ));
    end
% 
%     function grad = Grad(q)
%         C.Vertices = q;
%         grad = 2*Scal_Grad(C,C)- 2*Scal_Grad(C,S);
%     end

tgt.Grad = @Grad;


        

    function h = Plot(q,tag)
        if nargin<2
            tag = 'x-';
        end
        C.Vertices = q;
        h{1} = plotcurve(C,'g');
        h{2} = plotcurve(S,'r');
    end
tgt.Plot = @Plot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         gradient test       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function testgrad = TestGrad()
        
        J =  scal(C,C)+ scal(S,S) -2*scal(C,S);
        grad = Grad(C.Vertices);
        n = size(C.Vertices,1);
        for j = 1:10
            
            Y = rand(n,3);
            ScalGrad = sum(sum(grad.*Y));
            tab = [];
            
            for k = 1:10
                
                eps = 10^(-k);
                Ceps = C;
                Ceps.Vertices = Ceps.Vertices + eps*Y;
                Jeps =  scal(Ceps,Ceps)+scal(S,S)-2*scal(Ceps,S);
                Diff = (Jeps-J)/eps;
                tab = [tab,Diff - ScalGrad];
                
                
            end
            
            figure(j)
            clf
            plot(1:10,abs(tab))
        end
    end
tgt.testgrad = @TestGrad;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Some useful functions    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    

    
    function E = ComputeEdgesC3D(C)
        %compute the edges of all the triangles of the triangulation.
        E = C.Vertices(C.Faces(:,2),:)-C.Vertices(C.Faces(:,1),:);
    end

 
    function S = scal_prods(x,y)
        % computes matrix of <x(i,:),y(j,:)>
        [m,d] = size(x);
        [n,d] = size(y);
        S = zeros(m,n);
        for k = 1:d
            S = S + x(:,k)*y(:,k)';
        end
    end
    

    function totals = accummatrix(labels,X,sz)
        nCols = size(X,2);
        labels = [repmat(labels(:),nCols,1),kron(1:nCols,ones(1,numel(labels))).'];
        totals = accumarray(labels,X(:),sz);
    end

    function grad = gradient_scal(q)
        
        C.Vertices = q;
        C1 = C;
        %Extraction des elements utiles
        V1 = C1.Vertices;
        F1 = C1.Faces;
        E1 = ComputeEdgesC3D(C1);
        N = size(V1,1);
        n = size(E1,1);
        kerPP = kerp(V1,V1);
        normE1 = sqrt(sum(E1.*E1,2));
        NE1 = E1./repmat(normE1,1,3);
        d1kerpPP = permute(d1kerp(V1,V1),[2,3,1]);
        
        A = zeros(N,3);
        
        for d = 1:3
            A(:,d) = accumarray(F1(:,2),NE1(:,d),[N,1]) - accumarray(F1(:,1),NE1(:,d),[N,1]);
        end
        AA = scal_prods(A,A);
        grad1PP = sum(d1kerpPP.*repmat(AA,1,1,3),2);
        
        Atemp = permute(repmat(A,1,1,n),[3,1,2]);
        normE1tempA = repmat(normE1,1,N,3);
        NE1_Atemp = repmat(scal_prods(NE1,A),1,1,3);
        NE1tempA = permute(repmat(NE1,1,1,N),[1,3,2]);
        
        grad_edges_PP = (1./normE1tempA).*(Atemp - NE1_Atemp.*NE1tempA);
        grad_vertices_PP = zeros(N,N,3);
        
        for d = 1:3       
            
            KerPPF12 = kerPP(F1(:,2),:);
            KerPPF11 = kerPP(F1(:,1),:);
            
            grad_edges_kerPP2 = KerPPF12.*grad_edges_PP(:,:,d);
            grad_edges_kerPP1 = KerPPF11.*grad_edges_PP(:,:,d);
            
            grad_vertices_PP(:,:,d) = accummatrix(F1(:,2),grad_edges_kerPP2,[N,N]) + accummatrix(F1(:,1),grad_edges_kerPP1,[N,N]);
            grad_vertices_PP(:,:,d) = grad_vertices_PP(:,:,d) - accummatrix(F1(:,2),grad_edges_kerPP1,[N,N]) - accummatrix(F1(:,1),grad_edges_kerPP2,[N,N]);
            
        end
        
        grad2PP = sum(grad_vertices_PP,2);
        gradPP = grad1PP + grad2PP;
        grad = sigmaW^2*CoeffRenorm*2*pi^2/4*(squeeze(gradPP));
    end

tgt.Grad_scal = @gradient_scal;

end