clear

config{1}.x = [0;.5;0];
config{1}.mom = [1;0;0];
label{1} = '1lmk';

config{2}.x = [[-.4;.375;0],[.4;.625;0]];
config{2}.mom = [[1;0;0],[-1;0;0]];
label{2} = '2lmk';

config{3}.x = [[-.4;.5;0],[.4;.5;0]];
config{3}.mom = [[1;0;0],[-1;0;0]];
label{3} = 'headon';

Nx = 5;
label{4} = '5lmk';
config{4}.x = [rand(1,Nx)-.5;rand(1,Nx);zeros(1,Nx)];
config{4}.mom = randn(3,Nx);

sigmaV = .25;

% Nx = 100;
% s.x = rand(3,Nx);
% s.mom = randn(3,Nx);

for k=1:length(config)
    clear s
    s = config{k};
    mombase = s.mom;
    for a=1%-1:1
        for addb=0%[0,.03]
            i = 0;
            for coefmom=linspace(0,1,100)
                s.mom = coefmom*mombase;
                s.sigmaV = sigmaV;
                s.CppKer.Type = 'Tri';
                s.CppKer.FunctionTilde = 'WeightedGaussian';
                c = 1/s.sigmaV^2;
                if a<0
                    s.CppKer.CoefTilde = a;
                    s.CppKer.FunctionOrtho = 'WeightedGaussian';
                    b = addb-a/(2*c); % curl free
                    s.CppKer.CoefOrtho = b;
                else
                    s.CppKer.CoefTilde = a;
                    s.CppKer.FunctionOrtho = 'SpecGaussian';
                    dim = 2;
                    mu = dim/2-1;
                    b = addb+a*(2*mu+1)/(2*c);
                    b = 1;
                    s.CppKer.CoefOrtho(1) = b;
                    s.CppKer.CoefOrtho(2) = -a;
                end
                s.T = 100;
                
                s = shootingCpp(s);
                
                s.showtraj = 1;
                s.typefloat = 'double';
                s.optim_verbosemode = 1;
                s.transmatrix = eye(3);
                s.transvector = ones(3,1);
                s.showgrid = 1;
                s.gridsize = 100;
                s.usefgt = 0;
                s.tau = 1/(s.T-1);
                s.normcoefV = ones(1,s.T);
                s.target = {};
                s.show = {};
                figure(1)
                clf
                affiche(s);
                plot3(squeeze(s.X(1,:,:))',squeeze(s.X(2,:,:))',squeeze(s.X(3,:,:))','LineWidth',3)
                %plot([-.6,.6,.6,-.6],[0,0,1,1],'.r')
                axis equal
                axis off
                view(2)
                
                i = i+1;
                num = num2str(1000+i);num=num(2:end);
                ext = 'jpeg';
                print(['-d',ext],'-r200',['ess_',num,'.',ext])
            end
            eval(['!cp ess_',num,'.jpeg fig_',label{k},'_a=',num2str(a),'_b=',num2str(b),'.jpg'])
            eval(['!/usr/local/bin/convert -shave 150x150 ess*.jpeg anim_',label{k},'_a=',num2str(a),'_b=',num2str(b),'.gif'])
            !rm ess*.jpeg
            
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

        end
    end
end


% s.show = {0};
% s.target{1}.vx = 1:size(s.x,2);
% %s.target{1}.vy = 1:size(s.target{1}.y,2);
% makewrl('ess.wrl',s);


