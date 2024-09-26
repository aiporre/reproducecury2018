function [ButtonH,ButtonG,ButtonA] = gui_demo_matching(handles)

% Landmark matching example
% choose two landmark configurations x, y in R^(dn)
% and optimize either 1) or 2) :
% 1) J(p0) = gammaR E(p0,q0) + |y-q1|^2
% with q0=x, E(p,q) is the energy (hamiltonian) E(p,q)=<p,K(q)p>
% and q1 is solution at t=1 of Hamilton's equations for landmarks
% with initial conditions (p0,q0)
% 2) J(p0,q0) = gammaR E(p0,q0) + |x-q0|^2 + |y-q1|^2

N = 1;

% choose configurations
switch handles.config
    case 1
        % two landmarks crossing
        d = 2;
        n = 2;
        x = [0,0;1,.1];
        y = [1,0;0,.1];
    case 2
        % random configurations in 3D
        d = 3;
        n = 5;
        x = rand(n,d);
        y = x + randn(n,d)/10;
    case 3
        % circle to ellipse
        d = 2;
        n = 30;
        th = linspace(0,2*pi,n+1)';
        th = th(1:end-1);
        x = [cos(th),sin(th)];
        y = [2*cos(th),sin(th)];
    case 4
        % circle to some curve
        d = 2;
        n = 30;
        th = linspace(0,2*pi,n+1)';
        th = th(1:end-1);
        x = [cos(th),sin(th)];
        mod = 1+sin(4*th)/4;
        y = [2*cos(th).*mod,sin(th).*mod];
    case 5
        % circle to rotated ellipse
        d = 2;
        n = 30;
        th = linspace(0,2*pi,n+1)';
        th = th(1:end-1);
        x = 2*[cos(th),sin(th)];
        y = [4*cos(pi/4)*cos(th)+sin(pi/4)*sin(th)+1,cos(pi/4)*sin(th)-4*sin(pi/4)*cos(th)];
    case 6
        % circle to curve with measure matching
        d = 2;
        n = 30;
        th = linspace(0,2*pi,n+1)';
        th = th(1:end-1);
        x = .7*[cos(th),sin(th)];
        ny = 25;
        th = linspace(0,2*pi,ny+1)';
        th = th(1:end-1);
        mod = 1+sin(4*th)/8;
        y = [2*cos(th).*mod,sin(th).*mod];
    case 7
        % example of diffeomorphic shooting with landmarks
        % initial points and momenta
        q0 = [0,0;1,1];
        p0 = 5*[1,0;-1,0];
        d = size(q0,2); % dimension
        n = size(q0,1); % number of points
    case 8
        % Template estimation for landmarks
        N = 3; % number of configurations
        n = 5; % number of points for each configuration
        d = 3; % dimension
        x = rand(n,d);
        y = cell(1,N);
        for k=1:N
            y{k} = x + randn(n,d)/10;
        end
        q0 = rand(n,d);
    case 9
        d = 2;
        nstr = inputdlg('How many points','Custom configuration');
        n = str2num(nstr{1});
        x = [linspace(-1.5,1.5,n)',ones(n,1)];
        y = [linspace(-1.5,1.5,n)',-ones(n,1)];
        %         I = imread('ngc6543a.jpg');
        %         image(I,'XData',[-2,2],'YData',[-2,2])
        funGetx = dragpoints(x,'b');
        funGety = dragpoints(y,'r');
        axis([-2,2,-2,2])
        axis on
        ButtonD = uicontrol('Parent',gcf,'Style','pushbutton','String',...
            'Done','Units','normalized','Position',[0.85 0.03 0.1 0.05],...
            'Visible','on');
        set(ButtonD,'Callback',@DeleteMe)
        waitfor(msgbox('Move points (blue=source, red=target), then press "Done" button below','CreateMode','modal'))
        waitfor(ButtonD)
        x = funGetx();
        y = funGety();
    case 10
        % example of diffeomorphic shooting with landmarks
        d = 2;
        nstr = inputdlg('How many points','Custom configuration');
        n = str2num(nstr{1});
        q0 = (rand(n,2)-.5)*2.9;
        p0 = (rand(n,2)-.5);
        funGet = dragarrows(q0,p0,'r');
        axis([-2,2,-2,2])
        axis on
        ButtonD = uicontrol('Parent',gcf,'Style','pushbutton','String',...
            'Done','Units','normalized','Position',[0.85 0.03 0.1 0.05],...
            'Visible','on');
        set(ButtonD,'Callback',@DeleteMe)
        waitfor(msgbox('Move points and momentum vectors, then press "Done" button below','CreateMode','modal'))
        waitfor(ButtonD)
        [q0,p0] = funGet();
end

% choose kernel
sigma = handles.sigmaV; % kernel width
switch handles.kercode
    case 'matlab'
        ker = ScalarKernel(eval([handles.kerVtype,'Function']),sigma);
    case 'mex'
        ker = ScalarMexKernel(handles.kerVtype,sigma);
end

% kernel for measures
sigmaI = handles.sigmaI; % kernel width
switch handles.kercode
    case 'matlab'
        kerI = ScalarKernel(eval([handles.kerItype,'Function']),sigmaI);
    case 'mex'
        kerI = ScalarMexKernel(handles.kerItype,sigmaI);
end

% initialize positions and momenta
if handles.config<7 || handles.config==9
    q0 = x;
    p0{1} = zeros(n,d);
elseif handles.config == 8
    for k = 1:N
        p0{k} = zeros(n,d);
    end
elseif handles.config == 7
    p0 = p0 * sigma;
end

% gammaR parameter (0 for exact matching)
gammaR = handles.gamma;

% initialize matching
switch handles.typedef
    case 1
        for k = 1:N
            DM{k} = DiffeoMatch(HamiltSys(ker,n,d),gammaR);
        end
    case 2
        for k = 1:N
            DM{k} = SmallDefMatch(ker,gammaR);
        end
end
OD = OptimDiffeo(DM);

OD.SetOptimOption('Display','iter');
OD.SetOptimOption('HessUpdate','lbfgs');

if handles.config==8
    method = 2;
else
    method = 1;
end

if handles.config~=7 && handles.config~=10
    switch method
        case 1
            % optimize on initial momenta only
            if handles.config < 6 || handles.config==9
                DM{1}.AddTarget(TargetLandmarks(x,y));
            elseif handles.config == 6
                DM{1}.AddTarget(TargetMeasures(x,y,kerI));
            elseif handles.config == 8
                for k = 1:N
                    DM{k}.AddTarget(TargetLandmarks(x,y{k}));
                end
            end
            p0 = OD.Optimize_onlymom(p0,q0);
        case 2
            % optimize on initial positions and momenta
            if handles.config < 6
                DM.AddTarget(TargetLandmarks(x,x),0);
                DM.AddTarget(TargetLandmarks(x,y),1);
            elseif handles.config == 6
                DM.AddTarget(TargetMeasures(x,x,kerI),0);
                DM.AddTarget(TargetMeasures(x,y,kerI),1);
            elseif handles.config == 8
                for k = 1:N
                    DM{k}.AddTarget(TargetLandmarks(x,y{k}));
                end
            end
            [p0,q0] = OD.Optimize(p0,q0);
    end
end

% display results

[funAnim,funToggleGrid] = OD.Plot(p0,q0);
assignin('base','funAnim',funAnim);
assignin('base','funToggleGrid',funToggleGrid);
axis equal

ButtonH = uicontrol('Parent',gcf,'Style','pushbutton','String',...
    'View animation','Units','normalized','Position',[0.75 0.05 0.2 0.05],...
    'Visible','on','Callback','funAnim()');
ButtonG = uicontrol('Parent',gcf,'Style','checkbox','String',...
    'Show grid','Units','normalized','Position',[0.75 0.11 0.1 0.02],...
    'Value',1,'Visible','on','Callback','funToggleGrid()');
ButtonA = uicontrol('Parent',gcf,'Style','checkbox','String',...
    'Show axis','Units','normalized','Position',[0.85 0.11 0.1 0.02],...
    'Value',0,'Visible','on','Callback',@funToggleAxis);

function funToggleAxis(hObject, eventdata, handles)
if get(hObject,'Value')
    axis on
else
    axis off
end

function DeleteMe(hObject, eventdata, handles)
delete(hObject)
