function p = flowCpp(s,p,t1,t2)

if nargin<4
    t2 = s.T;
end

if nargin<3
    t1 = 1;
end

FileTemp1 = 'temp1.mch';
FileTemp2 = 'temp2.mch';
FileTemp3 = 'temp3.mch';

f = fopen(FileTemp1,'w');
fprintf(f,'Evol\n');
switch s.useDef
    case 'LargeDef'
        writeLargeDef(s,f);
    case 'LargeDef_InitParam'
        writeLargeDef_InitParam(s,f);
    case 'LargeDefSpec'
        writeLargeDefSpec(s,f);
    case 'SmallDef'
        writeSmallDef(s,f);
end
fclose(f);

f = fopen(FileTemp2,'w');
writeArr(p,f);
fclose(f);

pathtobins = [strrep(mfilename('fullpath'),mfilename,''),'../bin/'];
tagdim = num2str(size(s.x,1));

flowbin = [pathtobins,'flow',tagdim,'D'];

if ~exist(flowbin)
    disp(['recompiling c++ code for dimension ',tagdim])
    eval(['!../scriptcompile ',tagdim])
end
    
eval(['!',flowbin,' -f ',s.typefloat,' -m ',FileTemp1,' -d ',FileTemp2,' -o ',FileTemp3,' -b ',num2str(t1),' -e ',num2str(t2)])

f = fopen(FileTemp3,'r');
p = readArr(f);
fclose(f);
