
% essais pour TRI kernels with gaussian functions
% 2 cas envisagés :
% 1) example 1 du papier :
%   k_tilde(r) = a*exp(-c*r^2)  et  k_ortho(r) = (b-a*r^2)*exp(-c*r^2)
%   avec a >= 0 et b >= (d-1)*a/(2*c) avec d=2 ici
%   a=0 donne un noyau scalaire et b=(d-1)*a/(2*c) donne un noyau div-free
% 1) example 2 du papier, mais le coef a ici est l'opposé de celui du papier :
%   k_tilde(r) = a*exp(-c*r^2)  et  k_ortho(r) = b*exp(-c*r^2)
%   avec a <= 0 et b >= -a/(2*c)
%   a=0 donne un noyau scalaire et b=-a/(2*c) donne un noyau curl-free


clear

config{1}.x = [0;.5;0];
config{1}.mom = 10*[1;0;0];
config{1}.label = '1lmk';
config{1}.dec = [-0.5 -.1];
config{1}.scmom = .02;

config{2}.x = [[-.4;.375;0],[.4;.625;0]];
config{2}.mom = 80*[[1;0;0],[-1;0;0]];
config{2}.label = '2lmk';
config{2}.dec = [ -.75 -.35 ];
config{2}.scmom = .2;

config{3}.x = [[-.4;.5;0],[.4;.5;0]];
config{3}.mom = 10*[[1;0;0],[-1;0;0]];
config{3}.label = 'headon';
config{3}.dec = [ -.75 -.35 ];
config{3}.scmom = .2;

Nx = 5;
config{4}.label = '5lmk';
config{4}.x = [rand(1,Nx)-.5;rand(1,Nx);zeros(1,Nx)];
config{4}.mom = [7.5*randn(2,Nx);zeros(1,Nx)];
config{4}.dec = [ -.8 -.4 ];
config{4}.scmom = .2;

% list_config = config;
% clear config
% config{1} = list_config{2};

sigmaV = .25;

% Nx = 100;
% s.x = rand(3,Nx);
% s.mom = randn(3,Nx);

c = 1/sigmaV^2;
b = 1/(2*c);

%

for k=1:length(config)
    clear s
    s = config{k};
    for a=1%linspace(-1,1,5)
        i = 0;
        s.sigmaV = sigmaV;
        s.CppKer.Type = 'Tri';
        s.CppKer.FunctionTilde = 'WeightedGaussian';
        if a<0
            s.CppKer.CoefTilde = a;
            s.CppKer.FunctionOrtho = 'WeightedGaussian';
            s.CppKer.CoefOrtho = b;
        else
            s.CppKer.CoefTilde = a;
            s.CppKer.FunctionOrtho = 'SpecGaussian';
            s.CppKer.CoefOrtho(1) = b;
            s.CppKer.CoefOrtho(2) = -a;
        end
        s.T = 500;
        
        s = shootingCpp(s);
        
        s.showtraj = 1;
        s.typefloat = 'double';
        s.optim_verbosemode = 1;
        s.transmatrix = eye(3);
        s.transvector = ones(3,1);
        s.showgrid = 1;
        s.gridsize = 151;
        s.usefgt = 0;
        s.tau = 1/(s.T-1);
        s.normcoefV = ones(1,s.T);
        s.target = {};
        s.show = {};
        s.gridcolor = [1;0;0];
        
        figure(1)
        clf
        affiche(s);
        plot3(squeeze(s.X(1,:,1))',squeeze(s.X(2,:,1))',squeeze(s.X(3,:,1))','ko','MarkerSize',3)
        plot3(s.x(1,:),s.x(2,:),s.x(3,:),'ko','MarkerSize',3)
        %quiver3(s.x(1,:),s.x(2,:),s.x(3,:),s.mom(1,:,1),s.mom(2,:,1),s.mom(3,:,1),config{k}.scmom,'b','LineWidth',2)
        disp_mom = .01*s.mom(:,:,1);
        plot3(squeeze(s.X(1,:,:))',squeeze(s.X(2,:,:))',squeeze(s.X(3,:,:))','k','LineWidth',2)
        quiver3(s.x(1,:),s.x(2,:),s.x(3,:),disp_mom(1,:),disp_mom(2,:),disp_mom(3,:),0,'b','LineWidth',3)
        %plot([-.6,.6,.6,-.6],[0,0,1,1],'.r')
        axis equal
        axis off
        view(2)
        tag = ['a=',num2str(a),', b=',num2str(b),', c=',num2str(c)];
        switch a
            case -1
                tag = [tag,' (curl-free)'];
            case 0
                tag = [tag,' (scalar)'];
            case 1
                tag = [tag,' (div-free)'];
        end
        %plot(config{k}.dec(1)+[0,sigmaV],config{k}.dec(2)+[0,0],'k','LineWidth',7)
        rectangle('Position',[config{k}.dec(1),config{k}.dec(2),sigmaV,.01],...
            'LineWidth',1,'FaceColor',[0;1;0]) 
        title(tag)
        shg
        
        i = i+1;
        num = num2str(1000+i);num=num(2:end);
        ext = 'epsc';
        ext2 = 'eps';
        %print(['-d',ext],'-r200',['fig_',config{k}.label,'_a=',num2str(a),'_b=',num2str(b),'.',ext2])
        pause
    %eval(['!/usr/local/bin/convert -shave 150x150 ess*.jpeg anim_',label{k},'_a=',num2str(a),'_b=',num2str(b),'.gif'])
    %!rm ess*.jpeg
    
    %             % check area conservation
    %             n = 100;
    %             ech = linspace(0,1,n);
    %             p = [[-.5+ech,.5*ones(1,n),.5-ech,-.5*ones(1,n)];[zeros(1,n),ech,ones(1,n),1-ech]];
    %             p(3,1) = 0;
    %             pf = flowCpp(s,p);
    %             figure(2)
    %             clf
    %             hold on
    %             plot(p(1,:),p(2,:),'r')
    %             plot(pf(1,:),pf(2,:),'g')
    %             legend({['area inside=',num2str(polyarea(p(1,:),p(2,:)))],['area inside=',num2str(polyarea(pf(1,:),pf(2,:)))]})
     s.show = {0};
s.target{1}.vx = 1:size(s.x,2);
%s.target{1}.vy = 1:size(s.target{1}.y,2);
s.targetcolors = [0;0;0];
%makewrl(['fig_',config{k}.label,'_a=',num2str(a),'_b=',num2str(b),'.wrl'],s);

    end
   
end




