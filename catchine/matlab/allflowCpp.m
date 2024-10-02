function P = allflowCpp(s,p)

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
    case 'SmallDef'
        writeSmallDef(s,f);
end
fclose(f);

f = fopen(FileTemp2,'w');
writeArr(p,f);
fclose(f);

pathtobins = [strrep(mfilename('fullpath'),mfilename,''),'../bin/'];
tagdim = num2str(size(s.x,1));

allflowbin = [pathtobins,'allflow',tagdim,'D'];

if ~exist(allflowbin)
    disp(['recompiling c++ code for dimension ',tagdim])
    eval(['!../scriptcompile ',tagdim])
end
    
eval(['!',allflowbin,' -f ',s.typefloat,' -m ',FileTemp1,' -d ',FileTemp2,' -o ',FileTemp3])

f = fopen(FileTemp3,'r');
P = readArr(f);
fclose(f);