function writeSolver(f,s)

if ~isfield(s,'Solver')
        disp('using default solver : EulerTrapezoidal')
        s.Solver = 'EulerTrapezoidal';
end

fprintf(f,'%s\n',s.Solver);

