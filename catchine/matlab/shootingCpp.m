function s = shootingCpp(s)

% geodesic shooting for point matching

setenv('LD_LIBRARY_PATH',''); % for llull

FileTemp = 'temp.mch';
FileResults = 'resmatch.mch';

f = fopen(FileTemp,'w');
fprintf(f,'Evol\n');
writeLargeDef_InitParam(s,f);

fclose(f);

if ~isfield(s,'typefloat')
    disp('no typefloat field given, assuming type of float is double')
    s.typefloat = 'double';
end

pathtobins = [strrep(mfilename('fullpath'),mfilename,''),'../bin/'];
tagdim = num2str(size(s.x,1));

shootbin = [pathtobins,'shoot',tagdim,'D'];

if ~exist(shootbin)
    disp(['recompiling c++ code for dimension ',tagdim])
    eval(['!../scriptcompile ',tagdim])
end

eval(['!',shootbin,' ',s.typefloat,' ',FileTemp,' ',FileResults]);

typefloat = s.typefloat;
s = readResult(FileResults);
s.typefloat = typefloat;