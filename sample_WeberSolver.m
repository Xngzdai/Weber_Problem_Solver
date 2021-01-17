function sample_WeberSolver()
%   SAMPLE_WEBERSOLVER sample_WeberSolver() is a driver function for the
%   package Weber_Solver.m. It demonstrates how this package works and what
%   the input should be look like.
%
    W = [101 102 130 140 105];
    X = [19 92 93 94 59; 39 49 95 69 97];
    Y0 = [0;0];
    err = 0.01;
    
    loc = Weber_Solver(W,X,Y0,err);
    disp(loc)
end