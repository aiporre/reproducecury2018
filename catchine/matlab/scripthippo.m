

dosave = 0;

ExpLabel = 'Pt';
%ExpLabel = 'Hippo';

%KernelLabels = {'CurlFree','Scalar','DivFree'};
KernelLabels = {'Scalar'};

clear s target

for kk=1:length(KernelLabels)
    KernelLabel = KernelLabels{kk};
    s.T = 10;
    
    switch ExpLabel
        case 'Hippo'
    [s.x,target{1}.vx] = readbyu('s1002_hippoLsub500.byu');
    [target{1}.y,target{1}.vy] = readbyu('s1007_hippoLsub500.byu');
        case 'Pt'
    rf = .5;
    [V,F] = readbyu('001_lpt_rsub0_001_rpt_rsub0.byu'); [F,V] = reducepatch(F',V',rf); s.x = V'; target{1}.vx = F';
    %[V,F] = readbyu('003_lpt_rsub0_003_rpt_rsub0_proxmap.byu'); [F,V] = reducepatch(F',V',rf); s.x = V'; target{1}.vx = F';
    %S.Vertices = V; S.Faces = F; clf; surfplot(S,'r',1); camlight
    [V,F] = readbyu('007_lpt_rsub0_007_rpt_rsub0.byu'); [F,V] = reducepatch(F',V',rf); target{1}.y = V'; target{1}.vy = F';
    %[V,F] = readbyu('001_lpt_rsub0_001_rpt_rsub0.byu'); [F,V] = reducepatch(F',V',rf); target{1}.y = V'; target{1}.vy = F';
    %[V,F] = readbyu('003_lpt_rsub0_003_rpt_rsub0_proxmap.byu'); [F,V] = reducepatch(F',V',rf); target{1}.y = V'; target{1}.vy = F';
    %S.Vertices = V; S.Faces = F; hold on; surfplot(S,'g',1);
    end
    
    nfx = size(target{1}.vx,2);
    nfy = size(target{1}.vy,2);
    
    s.useDef = 'LargeDef';%_InitParam';    
        
    s.targetweights = 1;
    gammaR = 0;
    sigmaW = [5 1];
    optim_maxiter = [200,1000];
    optim_verbosemode = 2;
    s.optim_breakratio = 1e-6;
    s.optim_loopbreak = 10;
    s.rigidmatching = 0;
    
    
    target{1}.method = 'surfcurr';
    
    s.sigmaV = 10;
    
    switch computer
        case 'GLNXA64'
            target{1}.CppKer.Type = 'GaussGpu';
            s.CppKer.Type = 'GaussGpu';
        case 'MACI64'
            target{1}.CppKer.Type = 'SqDistScalar';
            target{1}.CppKer.Function = 'Gaussian';
            s.CppKer.Type = 'SqDistScalar';
            s.CppKer.Function = 'Gaussian';
    end
    
    % curl free kernel
    % s.CppKer = TRIkernel(KernelLabel,s.sigmaV);
    
    %target{1}.y = s.x+.2;%rand(3,3);
    %target{1}.vx = [1;2;3];
    %target{1}.vy = [1;2;3];
    target{1}.wx = ones(1,nfx);
    target{1}.wy = ones(1,nfy);
    
    nminims = length(sigmaW);
    for k=1:nminims
        target{1}.sigmaW = sigmaW(k);
        s.gammaR = gammaR;
        %s.useoptim = 'lbfgs_Quentin';
        s.optim_maxiter = optim_maxiter(k);
        s.optim_stepsize = optim_stepsize;
        s.optim_verbosemode = optim_verbosemode;
        s = matchCpp(s,target);
        
        %     clf;hold on;s.showpoints=0;s.show={'y','x','phi'};affiche(s);%h=plotsurf(Ssrc,'y');camlight;alpha(h,.5)
        %     pause
    end
    
    s.showpoints = 0;
    s.gridcolor = [0;0;0];
    s.gridmargin = .3;
    s.showgrid = 1;
    
    save([ExpLabel,KernelLabel,'.mat'])
    
    switch computer
        case 'MACI64'
            for k=0:3
                s.showtraj = ~~k;
                s.timeflow = k/3;
                clf
                h = affiche(s);
                axis equal
                axis off
                %view(88,62) % for hippos
                view(-9,-21) % for pts
                zoom(2)
                if dosave
                    saveas(gcf,[ExpLabel,KernelLabel,num2str(k),'o3.fig'])
                    print('-djpeg','-r400',[ExpLabel,KernelLabel,num2str(k),'o3.jpg'])
                    if k
                        delete(h.traj{1})
                        saveas(gcf,[ExpLabel,KernelLabel,num2str(k),'o3_notraj.fig'])
                        print('-djpeg','-r400',[ExpLabel,KernelLabel,num2str(k),'o3_notraj.jpg'])
                    end
                end
            end
    end
    
    % s.transmatrix = eye(3);
    % s.transvector = zeros(3,1);
    %
    %
    % %s.showgrid = 1;
    % makewrl('ess.wrl',s);
    
end



