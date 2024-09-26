function tgt = TargetNormalCycleCurves_Constant_linear(C,S,sigmaW,CoeffRenorm)

C.OrientedEdges = OrientedEdgesAdjacent(C);
S.OrientedEdges = OrientedEdgesAdjacent(S);

tgt.InitPos = C.Vertices;


%linear normal kernel on normal cycles for 3D curves

       

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
%Cylindric Scalar prod  %
%%%%%%%%%%%%%%%%%%%%%%%%%


function Scal = scal_constant(C1,C2)
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
        Scal = pi^2*sum(sum(kerPQ.*AB));
    end
    function scal = ScalVar3D(C1,C2,Phi)

    %compute the cylindric scalar product of two unions of segments   
        P1 = ComputeCenters(C1);
        P2 = ComputeCenters(C2);
        E1 = C1.Edges;
        E2 = C2.Edges;
        scal = sum(sum( kerp(P1,P2) .* scal_prods(E1,E2) .* cos(Phi) ));
    end


    function scal = Scal_sph_linear(C1,C2)
        V1 = C1.Vertices;
        V2 = C2.Vertices;
        N = size(V1,1);
        M = size(V2,1);
        T1 = C1.OrientedEdges;
        T2 = C2.OrientedEdges;
        nx = zeros(N,1);
        ny = zeros(1,M);
        for k = 1:N
            nx(k) = size(T1{k},1);
        end
        for k = 1:M
            ny(k) = size(T2{k},1);
        end
        scal = real(sum(sum(kerp(V1,V2).*(ones(N,M)-repmat(nx,1,M)/2).*(ones(N,M)-repmat(ny,N,1)/2))));
    end


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normal Cycle Scalar prod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Scal = scal(C1,C2)

        %liste des arretes normalisees
        C1.Edges = ComputeEdgesC3D(C1);
        C2.Edges = ComputeEdgesC3D(C2);
        C1.NormalizedEdges = Normalize(C1.Edges);
        C2.NormalizedEdges = Normalize(C2.Edges);

        %matrices des angles entre les arretes
        Phi = ComputeAngleEdges(C1.NormalizedEdges,C2.NormalizedEdges);
        Scal = CoeffRenorm*sigmaW^2*scal_constant(C1,C2) + pi^2/2*ScalVar3D(C1,C2,Phi) + 16*pi^2/3*CoeffRenorm*sigmaW^2*Scal_sph_linear(C1,C2);
    end
    scalS = scal(S,S);
    function A = Eval(q)
        C.Vertices = q;
        A = scal(C,C) + scalS -2*scal(C,S);
    end
tgt.Eval = @Eval;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Cylindric Gradkernel  %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function kernvargrad = KernelVarGrad()
        
    %fonction associee aux coefficients. Attention : il manque un produit associe aux aretes.   
        function ker = kern(phi)       
            ker = sin(phi); 
            inz = ker~=0;
            ker(inz) = ker(inz)./abs(sin(phi(inz)));
        end
        kernvargrad = @kern;
    end

kernvargrad = KernelVarGrad(); 

%%%%%%%%%%%%%%%%%%%%%%%%%
%    Cylindric  Grad    %
%%%%%%%%%%%%%%%%%%%%%%%%%

    function grad = Grad_constant(q)
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
        grad = 2*pi^2*(squeeze(gradPP)-squeeze(gradPQ));
    end

    function grad = GradVarC3D(C)
        
        %Extraction des elements utiles
        V1 = C.Vertices;
        E1 = C.Edges;
        E2 = S.Edges;
        P = ComputeCenters(C);
        Q = ComputeCenters(S);
        N = size(V1,1);
        n = size(E1,1);
        m = size(E2,1);
        N1 = Normalize(E1);
        N2 = Normalize(E2);
        %Matrices des angles entre les aretes
        Phi11 = ComputeAngleEdges(N1,N1);
        Phi12 = ComputeAngleEdges(N1,N2);
        T1 = OrientedEdgesAdjacent(C);
        %T2 = OrientedEdgesAdjacent(S);
        %remplissage des differentes matrices.
        
        %Pour la partie cylindrique, toutes les valeurs dependent soit du
        %milieu des aretes (pour la partie spatiale) soit des angles entre les
        %aretes (pour le noyau cylindrique). Dans les deux cas c'est une double
        %boucles sur les aretes.
        KernPP = cos(Phi11); %Matrice des valeurs du noyau cylindrique lors du produit scalaire de la courbe 1 avec elle meme
        KerpPP = kerp(P,P); %Matrice des valeurs du noyau spatial lors du produit scalaire de la courbe 1 avec elle meme
        d1kerpPP = d1kerp(P,P); %Matrice des gradients du noyau spatial lors du produit scalaire de la courbe 1 avec elle meme
        %d2kerpPP = d2kerp(P,P);
        KernvargradPP = kernvargrad(Phi11); %Matrice des gradients du noyau cylindrique lors du produit scalaire de la courbe 1 avec elle meme
        %On fait de meme pour le produit scalaire entre la courbe 1 et la
        %courbe 2.
        KernPQ = cos(Phi12);
        KerpPQ = kerp(P,Q);
        d1kerpPQ = d1kerp(P,Q);
        KernvargradPQ = kernvargrad(Phi12);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%Calcul du gradient%%%%%%%%%%%%%%%%%%%%%%%
        grad = zeros(N,3); %N lignes : nombre de sommets de la courbe 1
        %3 colones : valeur du gradient en chaque sommet.
        
        dim = size(d1kerpPQ,1);
        d1kerpPP = permute(d1kerpPP,[2,3,1]);
        d1kerpPQ = permute(d1kerpPQ,[2,3,1]);
        scE1 = scal_prods(E1,E1);
        scE12 = scal_prods(E1,E2);
        nE1 = sqrt(sum(E1.^2,2));
        oonE1 = (1./nE1)*ones(1,n);
        oonE1b = (1./nE1)*ones(1,m);
        oonE2 = (1./sqrt(sum(E2.^2,2)))*ones(1,n);
        coef1 = zeros(n,n,dim);
        coef2 = zeros(n,n,dim);
        coef3 = zeros(n,m,dim);
        coef4 = zeros(n,m,dim);
        
        coefa = scE1 .* KernPP;
        coefb = scE1.*KernvargradPP.*oonE1';
        coefc = (((1./nE1).^3)*ones(1,n)).*scE1;
        coefd = scE12 .* KernPQ;
        coefe = scE12.*KernvargradPQ.*oonE2';
        coeff = (((1./nE1).^3)*ones(1,m)).*scE12;
        for d = 1:dim
            E1d = E1(:,d)*ones(1,n);
            coef1(:,:,d) = d1kerpPP(:,:,d) .* coefa;        % + d2kerpPP{j,i})
            coef2(:,:,d) = 2 * KerpPP .* (KernPP .* E1d' + coefb.*(E1d'.*oonE1-coefc.*E1d));
            E2d = E2(:,d)*ones(1,n);
            E1d = E1(:,d)*ones(1,m);
            coef3(:,:,d) = d1kerpPQ(:,:,d) .* coefd;
            coef4(:,:,d) = 2 * KerpPQ .* (KernPQ .* E2d' + coefe.*(E2d'.*oonE1b-coeff.*E1d));
        end
        
        for k = 1:N %on construit le gradient pour chaque sommet k
            Ttemp = T1{k};%pour chaque sommet, on repertorie les aretes partant/arrivant sur ce sommet
            gradtemp = 0;
            for l = 1:size(Ttemp,1) %boucles sur les aretes de sommet k
                i = Ttemp(l)*sign(Ttemp(l));%arete non signee ayant k pour sommet
                gradtemp = gradtemp + sum(coef1(i,:,:) + coef2(i,:,:) * sign(Ttemp(l)),2) - sum(coef3(i,:,:) + coef4(i,:,:) * sign(Ttemp(l)),2);
            end
            %valeur du gradient pour le sommet k
            grad(k,:) = gradtemp;
        end
    end

    function grad = Grad_sph_linear(C)
        V1 = C.Vertices;
        V2 = S.Vertices;
        d1kerpPP = d1kerp(V1,V1);
        d1kerpPQ = d1kerp(V1,V2);
        d1kerpPP = permute(d1kerpPP,[2,3,1]);
        d1kerpPQ = permute(d1kerpPQ,[2,3,1]);
        N = size(V1,1);
        M = size(V2,1);
        T1 = C.OrientedEdges;
        T2 = S.OrientedEdges;
        nx = zeros(N,1);
        ny = zeros(1,M);
        for k = 1:N
            nx(k) = size(T1{k},1);
        end
        for k = 1:M
            ny(k) = size(T2{k},1);
        end
        A = (ones(N,N)-repmat(nx,1,N)/2).*(ones(N,N)-repmat(nx',N,1)/2);
        B = (ones(N,M)-repmat(nx,1,M)/2).*(ones(N,M)-repmat(ny,N,1)/2);
        
        DerivKerpPP = d1kerpPP.*repmat(A,1,1,3);
        DerivKerpPQ = d1kerpPQ.*repmat(B,1,1,3);
        
        grad = squeeze(2*(sum(DerivKerpPP,2)-sum(DerivKerpPQ,2)));

    end
%     function G = Grad(C)
%             d1kerp = d1gauss(sigmaW);
%             %Arete des courbes
%             C.Edges = ComputeEdgesC3D(C);
%             S.Edges = ComputeEdgesC3D(S);
%             G = 1/(sigmaW)^2*GradVarC3D(C);
%     end

    function G = Grad(q)
        G = zeros(size(q));
        C.Vertices = q;
        d1kerp = d1gauss(sigmaW);
        %Arete des courbes
        C.Edges = ComputeEdgesC3D(C);
        S.Edges = ComputeEdgesC3D(S);
        G = CoeffRenorm*sigmaW^2*Grad_constant(C.Vertices) + pi^2/2*GradVarC3D(C) + 16*pi^2/3*CoeffRenorm*sigmaW^2*Grad_sph_linear(C);
    end
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
% 
%     function h = Plot(tag)
%         if nargin<1
%             tag = 'x-';
%         end
%           h =  plot3([S.Vertices(S.Faces(:,1),1)';S.Vertices(S.Faces(:,2),1)'],[S.Vertices(S.Faces(:,1),2)';S.Vertices(S.Faces(:,2),2)'],[S.Vertices(S.Faces(:,1),3)';S.Vertices(S.Faces(:,2),3)'],'r','LineWidth',3);
%     end
% tgt.Plot = @Plot;

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


    function T = OrientedEdgesAdjacent(C)
    %pour chaque sommet k, T{k} liste les aretes dont k est un sommet, avec un
    %signe + si l orientation donnee au depart est rentrante, un - sinon.
        nvertices = size(C.Vertices,1);
        F = C.Faces;
        N = size(F,1);
        T = cell(nvertices,1);
        for i = 1:N
            T{F(i,1)} = [T{F(i,1)};-i];
            T{F(i,2)} = [T{F(i,2)};i];
        end
    end

    function N = Normalize(E)
    %Normalize the edges
        n = size(E,1);
        T = zeros(n,1);   
        for i = 1:n       
            T(i) = norm(E(i,:));   
        end    
        N = E./repmat(T,1,3);
    end

    function E = ComputeEdgesC3D(C)
        %compute the edges of all the triangles of the triangulation.
        E = C.Vertices(C.Faces(:,2),:)-C.Vertices(C.Faces(:,1),:);
    end

    function Phi = ComputeAngleEdges(NC,NS)
    %Matrice of all the angle betwenn two sets of normalized edges
        Phi = acos(NC*NS');
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

%     function L = ComputeLength(C1,C2)
%         E1 = C1.Edges;
%         E2 = C2.Edges;
%         L1 = sqrt(sum(E1.*E1,2));
%         L2 = sqrt(sum(E2.*E2,2));
%         L = L1*L2';
%     end

    function mu = ComputeCenters(S)
        %Compute the center of edges.
        mu = 1/2*(S.Vertices(S.Faces(:,1),:) + S.Vertices(S.Faces(:,2),:));
    end

    function totals = accummatrix(labels,X,sz)
        nCols = size(X,2);
        labels = [repmat(labels(:),nCols,1),kron(1:nCols,ones(1,numel(labels))).'];
        totals = accumarray(labels,X(:),sz);
    end
end

