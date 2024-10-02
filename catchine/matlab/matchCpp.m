function s = matchCpp(s,targets,tagdelay)

if nargin < 3 || ~strcmp(tagdelay,'delay')
    dodelay = 0;
else
    dodelay = 1;
end

if ~isfield(s,'typefloat')
    disp('no typefloat field given, assuming type of float is double')
    s.typefloat = 'double';
end


if ~isfield(s,'filename')
    s.filename = '';
end

FileTemp = [s.filename,'_temp.mch'];
FileResults = [s.filename,'_resmatch.mch'];

elapsedTime = -1;

if ~dodelay || ~exist(FileResults,'file')
    %    setenv('PATH',['/usr/local/cuda/bin:',getenv('PATH')]);
    %    setenv('DYLD_LIBRARY_PATH',['/usr/local/cuda/lib64:',getenv('DYLD_LIBRARY_PATH')]);
    %    setenv('LD_LIBRARY_PATH',[getenv('LD_LIBRARY_PATH') ':/usr/local/cuda/lib64']);
        setenv('LD_LIBRARY_PATH',''); % for pascal
    
    
    if isfield(s,'rigidmatching') && s.rigidmatching==1
        s.rigidmatching = 0;
        disp('Rigid matching not implemented in C++ code')
    end
    
    
    f = fopen(FileTemp,'w');
    fprintf(f,'Evol\n');
    if ~isfield(s,'useDef')
        s.useDef = 'LargeDef';
        disp('No useDef field given; assuming large deformations')
    end
    if ~isfield(s,'CppKer')
        if ~isfield(s,'usefgt') || s.usefgt==0
            s.CppKer.Type = 'SqDistScalar';
        else
            s.CppKer.Type = 'FastGauss';
            s.CppKer.order = 7;
            s.CppKer.K = 55;
            s.CppKer.e = 9;
        end
        s.CppKer.Function = 'Gaussian';
        disp('no CppKer field given (deformation). Assuming gaussian kernel')
    end
    switch s.useDef
        case 'LargeDef'
            writeLargeDef(s,f);
        case 'LargeDef_InitParam'
            writeLargeDef_InitParam(s,f);
        case 'LargeDefSpec'
            writeLargeDefSpec(s,f);
        case 'SmallDef'
            writeSmallDef(s,f);
        case 'FreeEvol'
            writeFreeEvol(s,f);
        otherwise
            error('unknown deformation type')
    end
    
    %    if (strcmp(s.CppKer.Type,'CauchyGpu') || strcmp(s.CppKer.Type,'GaussGpu')) && ~strcmp(s.typefloat,'float')
    %        s.typefloat = 'float';
    %        disp('redefining typefloat to ''float'' for use of Gpu code')
    %    end
    
    
    % writing the targets
    fprintf(f,'Targets\n');
    ntargets = length(targets);
    fprintf(f,'%d\n',ntargets);
    for i=1:ntargets
        if ~isfield(targets{i},'weight')
            targets{i}.weight = 1;
        end
        switch(targets{i}.method)
            case 'measures'
                if ~isfield(targets{i},'CppKer')
                    targets{i}.CppKer.Type = 'SqDistScalar';
                    targets{i}.CppKer.Function = 'Gaussian';
                    disp(['no CppKer field given for target ',num2str(i),'. Assuming scalar gaussian kernel'])
                end
                if isfield(targets{i},'sigmaI')
                     targets{i}.CppKer.Sigma = targets{i}.sigmaI;
                end
                fprintf(f,'Measure\n');
                fprintf(f,'Range=\n%d %d\n',targets{i}.vx(1),targets{i}.vx(end));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                fprintf(f,'Kernel=\n');
                writeKernel(f,targets{i}.CppKer)
                %                if (strcmp(targets{i}.CppKer.Type,'CauchyGpu')||strcmp(targets{i}.CppKer.Type,'GaussGpu')) && ~strcmp(s.typefloat,'float')
                %                    s.typefloat = 'float';
                %                    disp('redefining typefloat to ''float'' for use of Gpu code')
                %                end
                fprintf(f,'Y=\n');
                writeArr(targets{i}.y,f);
                fprintf(f,'WX,WY=\n');
                if ~isfield(targets{i},'wx')
                    nx = length(targets{i}.vx);
                    targets{i}.wx = ones(1,nx)/nx;
                end
                writeArr(targets{i}.wx(:)',f);
                if ~isfield(targets{i},'wy')
                    ny = size(targets{i}.y,2);
                    targets{i}.wy = ones(1,ny)/ny;
                end
                writeArr(targets{i}.wy(:)',f);
            case 'surfcurr'
                fprintf(f,'SurfCurr\n');
                fprintf(f,'Range=\n%d %d\n',min(targets{i}.vx(:)),max(targets{i}.vx(:)));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                if ~isfield(targets{i},'CppKer')
                    targets{i}.CppKer.Type = 'SqDistScalar';
                    targets{i}.CppKer.Function = 'Gaussian';
                    disp(['no CppKer field given for target ',num2str(i),'. Assuming scalar gaussian kernel'])
                end
                fprintf(f,'Kernel=\n');
                if isfield(targets{i},'sigmaW')
                    targets{i}.CppKer.Sigma = targets{i}.sigmaW;
                end
                writeKernel(f,targets{i}.CppKer)
                %               if (strcmp(targets{i}.CppKer.Type,'CauchyGpu')||strcmp(targets{i}.CppKer.Type,'GaussGpu')) && ~strcmp(s.typefloat,'float')
                %                   s.typefloat = 'float';
                %                   disp('redefining typefloat to ''float'' for use of Gpu code')
                %               end
                fprintf(f,'VY=\n');
                writeArr(targets{i}.y,f);
                fprintf(f,'FX=\n');
                writeArrInt(targets{i}.vx-min(targets{i}.vx(:))+1,f);
                fprintf(f,'FY=\n');
                writeArrInt(targets{i}.vy,f);
                fprintf(f,'WX,WY=\n');
                if ~isfield(targets{i},'wx')
                    targets{i}.wx = ones(1,size(targets{i}.vx,2));
                end
                writeArr(targets{i}.wx,f);
                if ~isfield(targets{i},'wy')
                    targets{i}.wy = ones(1,size(targets{i}.vy,2));
                end
                writeArr(targets{i}.wy,f);
                fprintf(f,'MatchingPursuitEpsilon=\n');
                if ~isfield(targets{i},'MatchingPursuitEpsilon')
                    targets{i}.MatchingPursuitEpsilon = 0;
                end
                fprintf(f,'%f\n',targets{i}.MatchingPursuitEpsilon);
            case 'curvecurr'
                fprintf(f,'CurveCurr\n');
                fprintf(f,'Range=\n%d %d\n',min(targets{i}.vx(:)),max(targets{i}.vx(:)));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                if ~isfield(targets{i},'CppKer.Type')
                    targets{i}.CppKer.Type = 'SqDistScalar';
                    targets{i}.CppKer.Function = 'Gaussian';
                    disp(['no CppKer field given for target ',num2str(i),'. Assuming scalar gaussian kernel'])
                end
                if isfield(targets{i},'sigmaW')
                    targets{i}.CppKer.Sigma = targets{i}.sigmaW;
                end
                fprintf(f,'Kernel=\n');
                writeKernel(f,targets{i}.CppKer)
                %               if (strcmp(targets{i}.CppKer.Type,'CauchyGpu')||strcmp(targets{i}.CppKer.Type,'GaussGpu')) && ~strcmp(s.typefloat,'float')
                %                  s.typefloat = 'float';
                %                   disp('redefining typefloat to ''float'' for use of Gpu code')
                %               end
                fprintf(f,'VY=\n');
                writeArr(targets{i}.y,f);
                fprintf(f,'FX=\n');
                writeArrInt(targets{i}.vx-min(targets{i}.vx(:))+1,f);
                fprintf(f,'FY=\n');
                writeArrInt(targets{i}.vy,f);
                fprintf(f,'WX,WY=\n');
                if ~isfield(targets{i},'wx')
                    targets{i}.wx = ones(1,size(targets{i}.vx,2));
                end
                writeArr(targets{i}.wx,f);
                if ~isfield(targets{i},'wy')
                    targets{i}.wy = ones(1,size(targets{i}.vy,2));
                end
                writeArr(targets{i}.wy,f);
                fprintf(f,'MatchingPursuitEpsilon=\n');
                if ~isfield(targets{i},'MatchingPursuitEpsilon')
                    targets{i}.MatchingPursuitEpsilon = 0;
                end
                fprintf(f,'%f\n',targets{i}.MatchingPursuitEpsilon);
            case 'curvecycle'
                fprintf(f,'CurveCycle\n');
                fprintf(f,'Range=\n%d %d\n',min(targets{i}.vx(:)),max(targets{i}.vx(:)));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                if ~isfield(targets{i},'CppKer')
                    targets{i}.CppKer.Type = 'SqDistScalar';
                    targets{i}.CppKer.Function = 'Gaussian';
                    disp(['no CppKer field given for target ',num2str(i),'. Assuming scalar gaussian kernel'])
                end
                if isfield(targets{i},'sigmaW')
                    targets{i}.CppKer.Sigma = targets{i}.sigmaW;
                end
                fprintf(f,'FunctionPoint=\n');
                fprintf(f,'%s,sigma=\n%f\n',targets{i}.CppKer.Function,targets{i}.CppKer.Sigma);
                %                if (strcmp(targets{i}.CppKer.Type,'CauchyGpu')||strcmp(targets{i}.CppKer.Type,'GaussGpu')) && ~strcmp(s.typefloat,'float')
                %                    s.typefloat = 'float';
                %                    disp('redefining typefloat to ''float'' for use of Gpu code')
                %                end
                if ~isfield(targets{i},'OrderConeFunction')
                    targets{i}.OrderConeFunction = 2;
                end
                fprintf(f,'OrderConeFunction=\n%f\n',targets{i}.OrderConeFunction);
                if ~isfield(targets{i},'AccConeFunction')
                    targets{i}.AccConeFunction = 1e-6;
                end
                fprintf(f,'AccConeFunction=\n%f\n',targets{i}.AccConeFunction);
                if ~isfield(targets{i},'OrderEdgeFunction')
                    targets{i}.OrderEdgeFunction = 2;
                end
                fprintf(f,'OrderEdgeFunction=\n%f\n',targets{i}.OrderEdgeFunction);
                if ~isfield(targets{i},'AccEdgeFunction')
                    targets{i}.AccEdgeFunction = 1e-6;
                end
                fprintf(f,'AccEdgeFunction=\n%f\n',targets{i}.AccEdgeFunction);
                fprintf(f,'VY=\n');
                writeArr(targets{i}.y,f);
                fprintf(f,'FX=\n');
                writeArrInt(targets{i}.vx-min(targets{i}.vx(:))+1,f);
                fprintf(f,'FY=\n');
                writeArrInt(targets{i}.vy,f);
            case 'curvacc'
                fprintf(f,'CurveAcc\n');
                fprintf(f,'Range=\n%d %d\n',min(targets{i}.vx(:)),max(targets{i}.vx(:)));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                if ~isfield(targets{i},'CppKer')
                    targets{i}.CppKer.Type = 'SqDistScalar';
                    targets{i}.CppKer.Function = 'Gaussian';
                    disp(['no CppKer field given for target ',num2str(i),'. Assuming scalar gaussian kernel'])
                end
                fprintf(f,'Kernel=\n');
                writeKernel(f,targets{i}.CppKer,targets{i}.sigmaW)
                %                if (strcmp(targets{i}.CppKer.Type,'CauchyGpu')||strcmp(targets{i}.CppKer.Type,'GaussGpu')) && ~strcmp(s.typefloat,'float')
                %                    s.typefloat = 'float';
                %                    disp('redefining typefloat to ''float'' for use of Gpu code')
                %                end
                fprintf(f,'VY=\n');
                writeArr(targets{i}.y,f);
                fprintf(f,'FX=\n');
                writeArrInt(targets{i}.vx-min(targets{i}.vx(:))+1,f);
                fprintf(f,'FY=\n');
                writeArrInt(targets{i}.vy,f);
                fprintf(f,'WX,WY=\n');
                if ~isfield(targets{i},'wx')
                    targets{i}.wx = ones(1,size(targets{i}.vx,2));
                end
                writeArr(targets{i}.wx,f);
                if ~isfield(targets{i},'wy')
                    targets{i}.wy = ones(1,size(targets{i}.vy,2));
                end
                writeArr(targets{i}.wy,f);
            case 'landmarks'
                fprintf(f,'Landmarks\n');
                fprintf(f,'Range=\n%d %d\n',targets{i}.vx(1),targets{i}.vx(end));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                fprintf(f,'Y=\n');
                writeArr(targets{i}.y,f);
            case 'landmarks_Norm1'
                fprintf(f,'Landmarks_Norm1\n');
                fprintf(f,'Range=\n%d %d\n',targets{i}.vx(1),targets{i}.vx(end));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                fprintf(f,'Y=\n');
                writeArr(targets{i}.y,f);
            case 'l2image'
                fprintf(f,'L2Image\n');
                fprintf(f,'Range=\n%d %d\n',targets{i}.rx(1),targets{i}.rx(end));
                fprintf(f,'Weight=\n%f\n',targets{i}.weight);
                fprintf(f,'SourceImage=\n');
                writeArr(targets{i}.imsource,f);
                fprintf(f,'TargetImage=\n');
                writeArr(targets{i}.imtarget,f);
                fprintf(f,'TargetGridBase=\n%f %f %f\n',targets{i}.basetarget);
                fprintf(f,'TargetGridVoxSize=\n%f %f %f\n',targets{i}.voxsizetarget);
        end
    end
    
    % Writing the Optimization method and parameters
    fprintf(f,'Optim\n');
    
    if isfield(s,'optim_useoptim')
        s.useoptim = s.optim_useoptim;
    end
    
    if ~isfield(s,'useoptim')
        s.useoptim = 'adaptdesc';
    end
    
    switch s.useoptim
        case 'fixedesc'
            WriteFixedDesc(s,f)
        case 'adaptdesc'
            WriteAdaptDesc(s,f)
        case 'lbfgs'
            WriteLbfgs_liblbfgs(s,f)
        case 'lbfgs_Quentin'
            WriteLbfgs_liblbfgs_Quentin(s,f)
        case 'lbfgs_dlib'
            WriteLbfgs_Dlib(s,f)
        case 'bfgs_dlib'
            WriteBfgs_Dlib(s,f)
        case 'cg_dlib'
            WriteCg_Dlib(s,f)
        otherwise
            error(['Unknown optimization method : ',s.useoptim]);
    end
    
    fclose(f);
    
    % type of float used in C++ code
    
    if ~isfield(s,'typefloat')
        s.typefloat = 'double';
    end
    
    pathtobins = [strrep(mfilename('fullpath'),mfilename,''),'../bin/'];
    tagdim = num2str(size(s.x,1));
    
    matchbin = [pathtobins,'match',tagdim,'D'];
        
    if ~exist(matchbin)
        disp(['recompiling c++ code for dimension ',tagdim])
        eval(['!../scriptcompile ',tagdim])
    end
    
    matchcmd = [matchbin,' -f ',s.typefloat,' -d ',FileTemp,' -o ',FileResults,' -v ',num2str(s.optim_verbosemode),...
        ' -w ',num2str(s.gammaR)];
    
    if dodelay
        s.matchcmd = matchcmd;
        s.delayed = 1;
    else
        eval(['!rm -f ',FileResults])
        tic;
        %eval(['!valgrind --leak-check=full ',command])
        eval(['!',matchcmd])
        %command
        %pause
        elapsedTime = toc;
    end
end
if ~dodelay || exist(FileResults,'file')
    typefloat = s.typefloat;
    s = readResult(FileResults);
    s.typefloat = typefloat;
    s.elapsedTime = elapsedTime;
    s.delayed = 0;
end



