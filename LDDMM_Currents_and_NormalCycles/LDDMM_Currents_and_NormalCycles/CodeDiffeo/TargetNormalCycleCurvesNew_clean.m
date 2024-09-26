function tgt = TargetNormalCycleCurvesNew_clean(C,S,sigmaW,TabCoeff,CoefRenormSph)

[nC,dC] = size(C.x);
if dC == 2
    C.x = [C.x,zeros(nC,1)];
    S.x = [S.x,zeros(size(S.x,1),1)];
end
C.Vertices = C.x;
C.Faces = C.G;
S.Vertices = S.x;
S.Faces = S.G;


lmax = length(TabCoeff)-1;
if nargin < 5
    CoefRenormSph = 1;
end

C.OrientedEdges = OrientedEdgesAdjacent(C);
S.OrientedEdges = OrientedEdgesAdjacent(S);
tgt.InitPos = C.Vertices;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EVALUATION OF THE NORMAL KERNEL k_n                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function eval = evalkernel(TabCoeff)
        function ker = ker(theta)
            ker = 0;
            for l = 0:lmax
                ker = ker + TabCoeff(l+1)*legendre(l,cos(theta));
            end
        end
        eval = @ker;
    end
tgt.EvalKernel = @evalkernel;
            
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CONSTRUCTION OF THE KERNELS (CYLINDRICAL AND SPHERICAL)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Cylindric Kernel    %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function kerncyl = KernelCyl(TabCoeff)
        %Calcul le noyau cylindrique. Depend uniquement de l angle phi entre deux
        %aretes.
        
        %Calcul des coefficients du noyau
        m = ceil(lmax/2);
        a = zeros(1,m);
        atemp = 0;
        
        for l = 0 : lmax
            atemp = atemp + TabCoeff(l+1) * (4 * ...
                IntLegendreSinCos(l,0,0,1)^2 + 2 * ...
                IntLegendreSinCos(l,2,0,1)^2 + 4 * ...
                IntLegendreSinCos(l,1,1,0)^2);
        end
        a(1) = atemp;
        if m > 1
        for k = 2:m
            atemp = 0;
            for l = 2*k-1:1:lmax
                atemp = atemp + TabCoeff(l+1) *( 4 * ...
                    IntLegendreSinCos(l,2*k-1,1,0)^2 + ...
                    2 * IntLegendreSinCos(l,2*k-2,0,1)^2 + ...
                    2 * IntLegendreSinCos(l,2*k,0,1)^2);
            end
            a(k) = atemp;
        end
        end
        %Calcul de la fonction de phi associee aux coefficients
        
        function ker = kern(phi)
            ker = zeros(size(phi));
            for i = 1:m
                ker = ker + a(i)*cos((2*i-1)*phi);
            end
        end
        kerncyl = @kern;
    end
    function CoeffCyl = coeffcyl(TabCoeff)
           m = ceil(lmax/2);
        a = zeros(1,m);
        atemp = 0;
        
        for l = 0 : lmax
            atemp = atemp + TabCoeff(l+1) * (4 * ...
                IntLegendreSinCos(l,0,0,1)^2 + 2 * ...
                IntLegendreSinCos(l,2,0,1)^2 + 4 * ...
                IntLegendreSinCos(l,1,1,0)^2);
        end
        a(1) = atemp;
        if m > 1
        for k = 2:m
            atemp = 0;
            for l = 2*k-1:1:lmax
                atemp = atemp + TabCoeff(l+1) *( 4 * ...
                    IntLegendreSinCos(l,2*k-1,1,0)^2 + ...
                    2 * IntLegendreSinCos(l,2*k-2,0,1)^2 + ...
                    2 * IntLegendreSinCos(l,2*k,0,1)^2);
            end
            a(k) = atemp;
        end
        end
        CoeffCyl = a;
    end
tgt.coeffcyl = @coeffcyl;
%%%%%%%%%%%%%%%%%%%%%%%%%
%Half Spherical Kernel  %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function kernHSHS = KernelHSHS(TabCoeff)
        %calcul du noyau demi spherique (le produit scalaire faisant apparaitre une
        %sphere complete est constant).
        
        %Coefficients du noyau
        m = ceil(lmax/2)+1;%+1 car le calcul du noyau spherique fait apparaitre un cos((m+1)phi)
        b = zeros(1,m);
        atemp = 0;
        TabCoeff = [TabCoeff,0,0];
        for l = 0 : lmax
            atemp = atemp + TabCoeff(l+1) * 4 * ...
                IntLegendreSinCos(l,0,2,0)^2 + TabCoeff(l+2) * 4 * ...
                IntLegendreSinCos(l+1,1,1,1)^2 + TabCoeff(l+3) * 2 * ...
                IntLegendreSinCos(l+2,2,2,0)^2;
        end
        b(1) = atemp;
        if m > 1
        for k = 2:m
            atemp = 0;
            for l = 2*k-1:1:lmax
                atemp = atemp + 1/((2*k-1)^2) * (TabCoeff(l+1) * 4 * ...
                    IntLegendreSinCos(l,2*k-1,1,1)^2 + TabCoeff(l) * ...
                    2 * IntLegendreSinCos(l-1,2*k-2,2,0)^2 + ...
                    TabCoeff(l+2) * 2 * IntLegendreSinCos(l+1,2*k,2,0)^2);
            end
            b(k) = atemp;
        end
        end
        %fonction de phi associee aux coefficients du noyau
        function ker = kern(phi)
            ker = 0;
            for i = 1:m
                ker = ker + b(i)*cos((2*i-1)*phi);
            end
        end
        kernHSHS = @kern;
        
    end

    function CoeffSph = coeffsph(TabCoeff)
        
       m = ceil(lmax/2)+1;%+1 car le calcul du noyau spherique fait apparaitre un cos((m+1)phi)
        b = zeros(1,m);
        atemp = 0;
        TabCoeff = [TabCoeff,0,0];
        for l = 0 : lmax
            atemp = atemp + TabCoeff(l+1) * 4 * ...
                IntLegendreSinCos(l,0,2,0)^2 + TabCoeff(l+2) * 4 * ...
                IntLegendreSinCos(l+1,1,1,1)^2 + TabCoeff(l+3) * 2 * ...
                IntLegendreSinCos(l+2,2,2,0)^2;
        end
        b(1) = atemp;
        if m > 1
        for k = 2:m
            atemp = 0;
            for l = 2*k-1:1:lmax
                atemp = atemp + 1/((2*k-1)^2) * (TabCoeff(l+1) * 4 * ...
                    IntLegendreSinCos(l,2*k-1,1,1)^2 + TabCoeff(l) * ...
                    2 * IntLegendreSinCos(l-1,2*k-2,2,0)^2 + ...
                    TabCoeff(l+2) * 2 * IntLegendreSinCos(l+1,2*k,2,0)^2);
            end
            b(k) = atemp;
        end
        end
        CoeffSph = b;
    end
tgt.coeffsph = @coeffsph;

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
kerncyl = KernelCyl(TabCoeff);
kernhs = KernelHSHS(TabCoeff);
d1kerp = d1gauss(sigmaW);
K = 0;
for l = 1:lmax+1
    K = K + pi^2*TabCoeff(l) *( 4 * IntLegendreSinCos(l-1,0,1,1)^2 + ...
         2 * IntLegendreSinCos(l-1,1,2,0)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Cylindric Scalar prod  %
%%%%%%%%%%%%%%%%%%%%%%%%%

    function scal = ScalNCC3Dcyl(C1,C2,Phi)
        
        %compute the cylindric scalar product of two unions of segments
        E1 = C1.Edges;
        E2 = C2.Edges;
        P1 = ComputeCenters(C1);
        P2 = ComputeCenters(C2);
        scal = sum(sum( kerp(P1,P2) .* scal_prods(E1,E2) .* kerncyl(Phi)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Scalar prod %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function scal = ScalNCC3Dsph(C1,C2,Phi)
        %calcul le produit scalaire spherique de deux unions de segments, C et S
        
        %K permet de calculer le produit scalaire ou une sphere entre jeu (que ce
        %soit avec une demi sphere ou avec une sphere
        %terme constant du noyau demi spherique
        azero = K/4;
        %noyau du produit scalaire entre les parties demi sph�riques
        V1 = C1.Vertices;
        V2 = C2.Vertices;
        N = size(V1,1);
        M = size(V2,1);
        T1 = C1.OrientedEdges;%OrientedEdgesAdjacent(C1);
        T2 = C2.OrientedEdges;%OrientedEdgesAdjacent(C2);
        F = C1.Faces;
        G = C2.Faces;
        n = size(F,1);
        m = size(G,1);
        A = zeros(N,M);
        a = kernhs(Phi);
        
        k1 = F(:,1)*ones(1,m);
        k2 = F(:,2)*ones(1,m);
        l1 = (G(:,1)*ones(1,n))';
        l2 = (G(:,2)*ones(1,n))';
        aplus = azero + a(:);
        amoins = azero - a(:);
        k1l1 = sub2ind(size(A),k1(:),l1(:));
        k1l2 = sub2ind(size(A),k1(:),l2(:));
        k2l1 = sub2ind(size(A),k2(:),l1(:));
        k2l2 = sub2ind(size(A),k2(:),l2(:));
        A2 = accumarray(k1l1,aplus,[M*N,1]);
        A2 = A2 + accumarray(k1l2,amoins,[M*N,1]);
        A2 = A2 + accumarray(k2l1,amoins,[M*N,1]);
        A2 = A2 + accumarray(k2l2,aplus,[M*N,1]);
        A = reshape(A2,N,M);
        
        nx = zeros(N,1);
        ny = zeros(1,M);
        for k = 1:N
            nx(k) = size(T1{k},1);
        end
        for k = 1:M
            ny(k) = size(T2{k},1);
        end
        scal = real(sum(sum(kerp(V1,V2).*((ones(N,M)-(repmat(nx,1,M)+repmat(ny,N,1))/2)*K+A))));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%
%Normal Cycle Scalar prod %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function Scal = scal(C1,C2)
        
        %liste des arretes normalisees
        C1.Edges = ComputeEdgesC3D(C1);
        C2.Edges = ComputeEdgesC3D(C2);
        C1.NormalizedEdges = Normalize(C1.Edges);
        C2.NormalizedEdges = Normalize(C2.Edges);
        %matrices des angles entre les arretes
        Phi = ComputeAngleEdges(C1.NormalizedEdges,C2.NormalizedEdges);
        Scal =  CoefRenormSph*sigmaW^2 * ScalNCC3Dsph(C1,C2,Phi) + ...
            ScalNCC3Dcyl(C1,C2,Phi);
    end

    function Scal = scaleval(q)
        C.Vertices = q;
	Scal = scal(C,S);
    end
tgt.Scal = @scaleval;


    function A = Eval(q)
        C.Vertices = q;
        A = scal(C,C) + scal(S,S) -2*scal(C,S);
    end
tgt.Eval = @Eval;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Cylindric Gradkernel  %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function kerncylgrad = KernelCylGrad(TabCoeff)
        %Calcul du gradient du noyau cylindrique
        
        %coefficient du noyau
        m = ceil(lmax/2);
        a = zeros(1,m);
        atemp = 0;
        
        for l = 0 : lmax
            atemp = atemp + TabCoeff(l+1) * (4 * ...
                IntLegendreSinCos(l,0,0,1)^2 + 2 * ...
                IntLegendreSinCos(l,2,0,1)^2 + 4 * ...
                IntLegendreSinCos(l,1,1,0)^2);
        end
        a(1) = atemp;
        if m > 1
        for k = 2:m
            atemp = 0;
            for l = 2*k-1:1:lmax
                atemp = atemp + TabCoeff(l+1) *( 4 * ...
                    IntLegendreSinCos(l,2*k-1,1,0)^2 + ...
                    2 * IntLegendreSinCos(l,2*k-2,0,1)^2 + ...
                    2 * IntLegendreSinCos(l,2*k,0,1)^2);
            end
            a(k) = atemp;
        end
        end
        %fonction associee aux coefficients. Attention : il manque un produit associe aux aretes.
        function ker = kern(phi)
            ker = zeros(size(phi));
            for i = 1:m
                ker = ker + (2*i-1)*a(i)*sin((2*i-1)*phi);
            end
            inz = ker~=0;
            ker(inz) = ker(inz)./abs(sin(phi(inz)));
        end
        kerncylgrad = @kern;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Gradkernel  %
%%%%%%%%%%%%%%%%%%%%%%%%%

    function kernHSHSgrad = KernelHSHSGrad(TabCoeff)
        %Calcul du gradient du noyau demi spherique
        
        %Coefficients des coefficients
        m = ceil(lmax/2)+1;%+1 car le calcul du noyau spherique fait apparaitre un cos((m+1)phi)
        b = zeros(1,m);
        atemp = 0;
        TabCoeff = [TabCoeff,0,0];
        for l = 0 : lmax
            atemp = atemp + TabCoeff(l+1) * 4 * ...
                IntLegendreSinCos(l,0,2,0)^2 + TabCoeff(l+2) * 4 * ...
                IntLegendreSinCos(l+1,1,1,1)^2 + TabCoeff(l+3) * 2 * ...
                IntLegendreSinCos(l+2,2,2,0)^2;
        end
        b(1) = atemp;
        if m > 1
        for k = 2:m
            atemp = 0;
            for l = 2*k-1:1:lmax
                atemp = atemp + 1/((2*k-1)^2) * (TabCoeff(l+1) * 4 * ...
                    IntLegendreSinCos(l,2*k-1,1,1)^2 + TabCoeff(l) * ...
                    2 * IntLegendreSinCos(l-1,2*k-2,2,0)^2 + ...
                    TabCoeff(l+2) * 2 * IntLegendreSinCos(l+1,2*k,2,0)^2);
            end
            b(k) = atemp;
        end
        end
        %function associee aux coefficients
        function ker = kern(phi)
            if phi==0
                ker = 0;
            else
                ker = 0;
                for i = 1:m
                    ker = ker + (2*i-1)*b(i)*sin((2*i-1)*phi);
                end
                inz = ker~=0;
                ker(inz) = ker(inz)./abs(sin(phi(inz)));
            end
        end
        kernHSHSgrad = @kern;
    end

kerncylgrad = KernelCylGrad(TabCoeff);
kernhsgrad = KernelHSHSGrad(TabCoeff);

%%%%%%%%%%%%%%%%%%%%%%%%%
%    Cylindric  Grad    %
%%%%%%%%%%%%%%%%%%%%%%%%%

    function grad = GradNCC3Dcyl(C)
        
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
        KerncylPP = kerncyl(Phi11); %Matrice des valeurs du noyau cylindrique lors du produit scalaire de la courbe 1 avec elle meme
        KerpPP = kerp(P,P); %Matrice des valeurs du noyau spatial lors du produit scalaire de la courbe 1 avec elle meme
        d1kerpPP = d1kerp(P,P); %Matrice des gradients du noyau spatial lors du produit scalaire de la courbe 1 avec elle meme
        %d2kerpPP = d2kerp(P,P);
        KerncylgradPP = kerncylgrad(Phi11); %Matrice des gradients du noyau cylindrique lors du produit scalaire de la courbe 1 avec elle meme
        %On fait de meme pour le produit scalaire entre la courbe 1 et la
        %courbe 2.
        KerncylPQ = kerncyl(Phi12);
        KerpPQ = kerp(P,Q);
        d1kerpPQ = d1kerp(P,Q);
        KerncylgradPQ = kerncylgrad(Phi12);
        
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
        
        coefa = scE1 .* KerncylPP;
        coefb = scE1.*KerncylgradPP.*oonE1';
        coefc = (((1./nE1).^3)*ones(1,n)).*scE1;
        coefd = scE12 .* KerncylPQ;
        coefe = scE12.*KerncylgradPQ.*oonE2';
        coeff = (((1./nE1).^3)*ones(1,m)).*scE12;
        for d = 1:dim
            E1d = E1(:,d)*ones(1,n);
            coef1(:,:,d) = d1kerpPP(:,:,d) .* coefa;        % + d2kerpPP{j,i})
            coef2(:,:,d) = 2 * KerpPP .* (KerncylPP .* E1d' + coefb.*(E1d'.*oonE1-coefc.*E1d));
            E2d = E2(:,d)*ones(1,n);
            E1d = E1(:,d)*ones(1,m);
            coef3(:,:,d) = d1kerpPQ(:,:,d) .* coefd;
            coef4(:,:,d) = 2 * KerpPQ .* (KerncylPQ .* E2d' + coefe.*(E2d'.*oonE1b-coeff.*E1d));
        end
        
        for k = 1:N %on construit le gradient pour chaque sommet k
            Ttemp = T1{k};%pour chaque sommet, on repertorie les aretes partant/arrivant sur ce sommet
            gradtemp = 0;
            for l = 1:size(Ttemp,1) %boucles sur les aretes de sommet k
                i = Ttemp(l)*sign(Ttemp(l));%arete non signee ayant k pour sommet
                gradtemp = gradtemp + sum(coef1(i,:,:) + coef2(i,:,:) * ...
                    sign(Ttemp(l)),2) - sum(coef3(i,:,:) + coef4(i,:,:) * ...
                    sign(Ttemp(l)),2);
            end
            %valeur du gradient pour le sommet k
            grad(k,:) = gradtemp;
        end
    end

   
%%%%%%%%%%%%%%%%%%%%%%%%%
%    Spherical Grad     %
%%%%%%%%%%%%%%%%%%%%%%%%%

    function grad = GradNCC3Dsph(C)
        %extraction des donnees utiles pour le calcul du gradient
        V1 = C.Vertices;
        dim = size(V1,2);
        V2 = S.Vertices;
        E1 = C.Edges;
        E2 = S.Edges;
        F1 = C.Faces;
        F2 = S.Faces;
        N = size(V1,1);
        M = size(V2,1);
        n = size(E1,1);
        m = size(E2,1);
        N1 = Normalize(E1);
        N2 = Normalize(E2);
        Phi11 = ComputeAngleEdges(N1,N1);
        Phi12 = ComputeAngleEdges(N1,N2);
        T1 = C.OrientedEdges; %Adjacent(C);
        T2 = S.OrientedEdges;%Adjacent(S);
        azero = K/4;
        KernhsPP = kernhs(Phi11);
        KernhsgradPP = kernhsgrad(Phi11);
        KerpPP = kerp(V1,V1);
        d1kerpPP = d1kerp(V1,V1);
        %d2kerpPP = d2kerp(V1,V1);
        KernhsPQ = kernhs(Phi12);
        KernhsgradPQ = kernhsgrad(Phi12);
        KerpPQ = kerp(V1,V2);
        d1kerpPQ = d1kerp(V1,V2);
        A = zeros(N,N);
        Agrad = zeros(N,N,dim);
        B = zeros(N,M);
        Bgrad = zeros(N,M,dim);
        d1kerpPP = permute(d1kerpPP,[2,3,1]);
        d1kerpPQ = permute(d1kerpPQ,[2,3,1]);
        scE1 = scal_prods(E1,E1);
        scE12 = scal_prods(E1,E2);
        nE1 = sqrt(sum(E1.^2,2));
        nE2 = sqrt(sum(E2.^2,2));
        oonE1 = (1./nE1)*ones(1,n);
        oonE1b = (1./nE1)*ones(1,m);
        oonE2 = (1./nE2)*ones(1,n);
        coef1 = zeros(n,n,dim);
        coef2 = zeros(n,m,dim);
        coefa = (((1./nE1).^3)*ones(1,n)).*scE1;
        coefb = (((1./nE1).^3)*ones(1,m)).*scE12;
        %gradient spherique de la norme de la difference des cycles normaux
        %associes aux courbes C1 et C2.
        grad = zeros(N,3);
        for d = 1:dim
            E1d = E1(:,d)*ones(1,n);
            E2d = E2(:,d)*ones(1,n);
            E1dm = E1(:,d)*ones(1,m);
            coef1(:,:,d) = KernhsgradPP.*oonE1'.*(oonE1.*E1d'-coefa.*E1d);
            coef2(:,:,d) = KernhsgradPQ.*oonE2'.*(oonE1b.*E2d'-coefb.*E1dm);
        end
        
        k1 = F1(:,1)*ones(1,n);
        k2 = F1(:,2)*ones(1,n);
        l1 = (F1(:,1)*ones(1,n))';
        l2 = (F1(:,2)*ones(1,n))';
        Kplus = azero + KernhsPP(:);
        Kmoins = azero - KernhsPP(:);
        k1l1 = sub2ind(size(A),k1(:),l1(:));
        k1l2 = sub2ind(size(A),k1(:),l2(:));
        k2l1 = sub2ind(size(A),k2(:),l1(:));
        k2l2 = sub2ind(size(A),k2(:),l2(:));
        At = accumarray(k1l1,Kplus,[N*N,1]);
        At = At + accumarray(k1l2,Kmoins,[N*N,1]);
        At = At + accumarray(k2l1,Kmoins,[N*N,1]);
        At = At + accumarray(k2l2,Kplus,[N*N,1]);
        A = reshape(At,N,N);
        for d = 1:dim
            tmp = (KerpPP(k1l1)-KerpPP(k2l1)).*reshape(coef1(:,:,d),n*n,1);
            Ag = accumarray(k1l1,-tmp,[N*N,1]);
            Ag = Ag + accumarray(k2l1,tmp,[N*N,1]);
            tmp = (KerpPP(k1l2)-KerpPP(k2l2)).*reshape(coef1(:,:,d),n*n,1);
            Ag = Ag + accumarray(k1l2,tmp,[N*N,1]);
            Ag = Ag + accumarray(k2l2,-tmp,[N*N,1]);
            Agrad(:,:,d) = reshape(Ag,N,N);
        end
        
        k1 = F1(:,1)*ones(1,m);
        k2 = F1(:,2)*ones(1,m);
        l1 = (F2(:,1)*ones(1,n))';
        l2 = (F2(:,2)*ones(1,n))';
        Kplus = azero + KernhsPQ(:);
        Kmoins = azero - KernhsPQ(:);
        k1l1 = sub2ind(size(B),k1(:),l1(:));
        k1l2 = sub2ind(size(B),k1(:),l2(:));
        k2l1 = sub2ind(size(B),k2(:),l1(:));
        k2l2 = sub2ind(size(B),k2(:),l2(:));
        Bt = accumarray(k1l1,Kplus,[N*M,1]);
        Bt = Bt + accumarray(k1l2,Kmoins,[N*M,1]);
        Bt = Bt + accumarray(k2l1,Kmoins,[N*M,1]);
        Bt = Bt + accumarray(k2l2,Kplus,[N*M,1]);
        B = reshape(Bt,N,M);
        for d = 1:dim
            tmp = (KerpPQ(k1l1)-KerpPQ(k2l1)).*reshape(coef2(:,:,d),n*m,1);
            Bg = accumarray(k1l1,-tmp,[N*M,1]);
            Bg = Bg + accumarray(k2l1,tmp,[N*M,1]);
            tmp = (KerpPQ(k1l2)-KerpPQ(k2l2)).*reshape(coef2(:,:,d),n*m,1);
            Bg = Bg + accumarray(k1l2,tmp,[N*M,1]);
            Bg = Bg + accumarray(k2l2,-tmp,[N*M,1]);
            Bgrad(:,:,d) = reshape(Bg,N,M);
        end
        nx = zeros(N,1);
        ny = zeros(1,M);
        
        for k = 1:N
            nx(k) = size(T1{k},1);
        end
        
        for l = 1:M
            ny(l) = size(T2{l},1);
        end
        A = (ones(N,N)-(repmat(nx,1,N)+repmat(nx',N,1))/2)*K+A;
        B = (ones(N,M)-(repmat(nx,1,M)+repmat(ny,N,1))/2)*K+B;
        DerivKerpPP = d1kerpPP.*repmat(A,1,1,3);
        DerivKerpPQ = d1kerpPQ.*repmat(B,1,1,3);
        for k = 1:N
            gradtemp =2*( sum(Agrad(k,:,:) + DerivKerpPP(k,:,:),2) - ...
                sum(DerivKerpPQ(k,:,:)+Bgrad(k,:,:),2 ));
            grad(k,:) = gradtemp;
        end
    end



    function G = Grad(q)
        G = zeros(size(q));
        C.Vertices = q;
        d1kerp = d1gauss(sigmaW);
        %Arete des courbes
        C.Edges = ComputeEdgesC3D(C);
        S.Edges = ComputeEdgesC3D(S);
        G = CoefRenormSph* sigmaW^2 * GradNCC3Dsph(C) + GradNCC3Dcyl(C);
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
            
            for k = 1:6
                
                eps = 10^(-k);
                Ceps = C;
                Ceps.Vertices = Ceps.Vertices + eps*Y;
                Jeps =  scal(Ceps,Ceps)+scal(S,S)-2*scal(Ceps,S);
                Diff = (Jeps-J)/eps;
                tab = [tab,Diff - ScalGrad];
                
                
            end
            
            figure(j)
            clf
            plot(1:6,abs(tab))
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

    function f = factorielstep(m,k)
        f = 1;
        while m > 0
            f = f*m;
            m = m-k;
        end
    end

    function I = IntSinCos(k,l)
        %compute int_0^pi sin^k theta cos^j theta
        if mod(l,2)== 1
            I = 0;
        else if (mod(k,2) == 0)
                I = factorielstep(l-1,2) * factorielstep(k-1,2) * ...
                    factorielstep(k+l,1) * pi/(factorielstep(k+l,2)^2 * ...
                    factorielstep(k+l-1,2));
            else I = 2 * factorielstep(l-1,2) * factorielstep(k-1,2) * ...
                    factorielstep(k+l-1,2) / (factorielstep(k+l,1));
            end
        end
    end

    function test = IntSinCosTest(k,l)
        f= @(x) sin(x).^k.*cos(x).^l;
        test = IntSinCos(k,l)-integral(f,0,pi)
    end
tgt.testsin = @IntSinCosTest;

    function I = IntLegendreSinCos(l,m,i,j)
        % Compute int_0^pi P_l^m(cos theta)sin^i theta cos^j theta
        k = ceil((l+m)/2);
        I = 0;
        if l < m
            I = 0;
        elseif (l ==0) && (m==0)
            I = (1/sqrt(4*pi))*IntSinCos(i,j);
        elseif (m == 0) && (l~= 0)
            while (k <= l)
                I = I + (-1)^(l-k) * factorielstep(2*k-1,2) / ...
                    (factorial(l-k) * factorial(2*k-l) * 2^(l-k)) * ...
                    IntSinCos(i,2*k-l+j);
                k = k+1;
            end
            I = real(sqrt((2*l+1)/(4*pi))*I);
        else
            while (k <= l)
                I = I + (-1)^(k+l+m) * factorielstep(2*k-1,2) / ...
                    (factorial(l-k) * factorial(2*k-(m+l)) * 2^(l-k)) * ...
                    IntSinCos(m+i,2*k-m-l+j);
                k = k+1;
            end
            I = real(sqrt((2*l+1)*factorial(l-m)/(factorial(l+m)*2*pi))*I);
        end
    end
tgt.intlegsincos = @IntLegendreSinCos;
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
    function mu = ComputeCenters(S)
        %Compute the center of edges.
        mu = 1/2*(S.Vertices(S.Faces(:,1),:) + S.Vertices(S.Faces(:,2),:));
    end
end
