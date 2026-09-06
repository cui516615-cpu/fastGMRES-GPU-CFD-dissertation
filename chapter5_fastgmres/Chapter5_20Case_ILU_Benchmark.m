clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf(' FINAL 20-CASE BENCHMARK FOR CHAPTER 5\n');
fprintf(' 10 MATRICES: WITHOUT + WITH ILU(0) PRECONDITIONING\n');
fprintf('============================================================\n\n');


%% ========================================================================
% 1. Matrix files
% ========================================================================

files = {
    'naca12A.mat'
    'ML_Laplace.mat'
    'vas_stokes_1M.mat'
    'Transport.mat'
    'cavity10.mat'
    'ns3Da.mat'
    'atmosmodl.mat'
    'atmosmodd.mat'
    'GT01R.mat'
    'poisson3Da.mat'
};

matrix_names = {
    'NACA12'
    'ML_Laplace'
    'vas_stokes_1M'
    'Transport'
    'cavity10'
    'ns3Da'
    'atmosmodl'
    'atmosmodd'
    'GT01R'
    'poisson3Da'
};

methods = {
    'GMRES(50)'
    'BiCGSTAB'
    'fastGMRES'
};

numMatrices = length(files);

% Each matrix is tested:
%
%   1. without preconditioning
%   2. with ILU(0) right preconditioning
%
numCases = 2*numMatrices;


%% ========================================================================
% 2. Check that all matrix files exist
% ========================================================================

fprintf('Checking matrix files...\n');

for p = 1:numMatrices

    if ~isfile(files{p})

        error( ...
            'Matrix file "%s" was not found in the current MATLAB folder.', ...
            files{p});

    end

    fprintf('  [%2d/10] %-22s : OK\n', ...
        p, files{p});

end

fprintf('\nAll 10 matrix files were found.\n\n');


%% ========================================================================
% 3. COMMON solver parameters
% ========================================================================

tol = 1e-7;


% -------------------------------------------------------------------------
% GMRES
% -------------------------------------------------------------------------

restart = 50;

% Maximum Krylov steps:
%
%       40 * 50 = 2000
%

gmres_max_outer = 40;


% -------------------------------------------------------------------------
% BiCGSTAB
% -------------------------------------------------------------------------

bicg_max_iter = 1000;


% -------------------------------------------------------------------------
% fastGMRES
% -------------------------------------------------------------------------

fast_max_outer = 150;

opts = struct;

opts.verbose = 0;

% Maximum wall-clock protection for one fastGMRES solve
opts.maxtime = 600;

% Pure non-orthogonal power basis
opts.t = 0;

% Maximum inner sGMRES dimension
opts.m = 80;

% Sketch dimension
opts.s = 2*(opts.m + 1);

% Conditioning safeguard
opts.cndtol = 1e15;


fprintf('============================================================\n');
fprintf(' COMMON SOLVER CONFIGURATION\n');
fprintf('============================================================\n');

fprintf('Target tolerance        = %.1e\n', ...
    tol);

fprintf('\nGMRES restart           = %d\n', ...
    restart);

fprintf('GMRES maximum steps     = %d\n', ...
    restart*gmres_max_outer);

fprintf('\nBiCGSTAB max iterations = %d\n', ...
    bicg_max_iter);

fprintf('\nfastGMRES outer maxit   = %d\n', ...
    fast_max_outer);

fprintf('fastGMRES inner m       = %d\n', ...
    opts.m);

fprintf('fastGMRES sketch s      = %d\n', ...
    opts.s);

fprintf('fastGMRES truncation t  = %d\n', ...
    opts.t);

fprintf('fastGMRES cndtol        = %.1e\n', ...
    opts.cndtol);

fprintf('fastGMRES maxtime       = %.0f s\n\n', ...
    opts.maxtime);


%% ========================================================================
% 4. ILU configuration
%
% Explicit zero-fill incomplete LU:
%
%       A approximately L*U
%
% This corresponds to ILU(0).
%
% IMPORTANT:
%
% The factorization time is recorded separately and is NOT included in
% the individual iterative solver runtimes.
% ========================================================================

ilu_setup = struct;

ilu_setup.type = 'nofill';

ilu_setup.milu = 'off';


fprintf('============================================================\n');
fprintf(' ILU PRECONDITIONER CONFIGURATION\n');
fprintf('============================================================\n');

fprintf('ILU type                = nofill (ILU(0))\n');

fprintf('Modified ILU            = off\n');

fprintf('Orientation              = right preconditioning\n');

fprintf('ILU setup time           = recorded separately\n\n');


%% ========================================================================
% 5. Check fastgmres.m
% ========================================================================

if exist('fastgmres', 'file') ~= 2

    error( ...
        'fastgmres.m was not found in the current MATLAB folder or path.');

end


%% ========================================================================
% 6. Allocate metadata arrays for all 20 cases
%
% For each matrix:
%
%       odd row  = no preconditioning
%       even row = ILU(0) preconditioning
%
% Example:
%
%       row 1 : NACA12
%       row 2 : NACA12-p
%
%       row 3 : ML_Laplace
%       row 4 : ML_Laplace-p
% ========================================================================

Case = strings(numCases,1);

Matrix = strings(numCases,1);

Mode = strings(numCases,1);

N = zeros(numCases,1);

NNZ = zeros(numCases,1);

RHS_Source = strings(numCases,1);

ILU_Status = strings(numCases,1);

ILU_Setup_Time = NaN(numCases,1);


%% ========================================================================
% 7. Allocate GMRES results
% ========================================================================

GMRES_Time = NaN(numCases,1);

GMRES_RelRes = Inf(numCases,1);

GMRES_Steps = NaN(numCases,1);

GMRES_Converged = false(numCases,1);


%% ========================================================================
% 8. Allocate BiCGSTAB results
% ========================================================================

BiCG_Time = NaN(numCases,1);

BiCG_RelRes = Inf(numCases,1);

BiCG_Iter = NaN(numCases,1);

BiCG_Converged = false(numCases,1);


%% ========================================================================
% 9. Allocate fastGMRES results
% ========================================================================

fast_Time = NaN(numCases,1);

fast_RelRes = Inf(numCases,1);

fast_OuterIter = NaN(numCases,1);

fast_Converged = false(numCases,1);


%% ========================================================================
% 10. Runtime matrix used for performance profiles
%
% Column 1 = GMRES(50)
% Column 2 = BiCGSTAB
% Column 3 = fastGMRES
%
% Inf means target tolerance was not achieved.
% ========================================================================

performance_time = Inf(numCases,3);


%% ========================================================================
% 11. Main benchmark loop
% ========================================================================

for p = 1:numMatrices

    fprintf('\n\n');
    fprintf('============================================================\n');
    fprintf(' MATRIX %d / %d : %s\n', ...
        p, numMatrices, matrix_names{p});
    fprintf('============================================================\n');


    %% --------------------------------------------------------------------
    % Row numbers for current matrix
    % ---------------------------------------------------------------------

    row_no_prec = 2*p - 1;

    row_ilu = 2*p;


    %% --------------------------------------------------------------------
    % Load matrix
    % ---------------------------------------------------------------------

    fprintf('\n[1] Loading matrix...\n');

    [A,b,rhs_source] = ...
        load_matrix_problem(files{p},1000+p);

    A = sparse(A);

    b = full(b(:));

    n = size(A,1);


    %% --------------------------------------------------------------------
    % Basic checks
    % ---------------------------------------------------------------------

    if size(A,1) ~= size(A,2)

        error('%s is not square.', ...
            files{p});

    end

    if length(b) ~= n

        error('RHS dimension mismatch for %s.', ...
            files{p});

    end

    if norm(b) == 0

        error('RHS is zero for %s.', ...
            files{p});

    end


    %% --------------------------------------------------------------------
    % Metadata for BOTH cases
    % ---------------------------------------------------------------------

    Case(row_no_prec) = ...
        string(matrix_names{p});

    Case(row_ilu) = ...
        string(matrix_names{p}) + "-p";


    Matrix(row_no_prec) = ...
        string(matrix_names{p});

    Matrix(row_ilu) = ...
        string(matrix_names{p});


    Mode(row_no_prec) = ...
        "No preconditioner";

    Mode(row_ilu) = ...
        "ILU(0) right";


    N(row_no_prec) = n;

    N(row_ilu) = n;


    NNZ(row_no_prec) = nnz(A);

    NNZ(row_ilu) = nnz(A);


    RHS_Source(row_no_prec) = rhs_source;

    RHS_Source(row_ilu) = rhs_source;


    ILU_Status(row_no_prec) = ...
        "N/A";


    fprintf('Dimension       : %d x %d\n', ...
        n,n);

    fprintf('NNZ             : %d\n', ...
        nnz(A));

    fprintf('Matrix density  : %.6e\n', ...
        nnz(A)/(double(n)*double(n)));

    fprintf('RHS             : %s\n', ...
        rhs_source);


    %% --------------------------------------------------------------------
    % Estimate fastGMRES main basis memory
    % ---------------------------------------------------------------------

    total_columns = ...
        (fast_max_outer+1) + ...
        fast_max_outer + ...
        (opts.m+1);

    basis_GB = ...
        8*double(n)*double(total_columns)/(1024^3);

    fprintf( ...
        'Estimated main fastGMRES basis memory: %.2f GB\n', ...
        basis_GB);


    %% ====================================================================
    % CASE A:
    % WITHOUT PRECONDITIONING
    % =====================================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' CASE A: %s -- WITHOUT PRECONDITIONING\n', ...
        matrix_names{p});
    fprintf('============================================================\n');


    x0 = zeros(n,1);


    R = run_solver_case( ...
        A, ...
        A, ...
        @(y) y, ...
        b, ...
        x0, ...
        tol, ...
        restart, ...
        gmres_max_outer, ...
        bicg_max_iter, ...
        fast_max_outer, ...
        opts);


    GMRES_Time(row_no_prec) = ...
        R.GMRES_Time;

    GMRES_RelRes(row_no_prec) = ...
        R.GMRES_RelRes;

    GMRES_Steps(row_no_prec) = ...
        R.GMRES_Steps;

    GMRES_Converged(row_no_prec) = ...
        R.GMRES_Converged;


    BiCG_Time(row_no_prec) = ...
        R.BiCG_Time;

    BiCG_RelRes(row_no_prec) = ...
        R.BiCG_RelRes;

    BiCG_Iter(row_no_prec) = ...
        R.BiCG_Iter;

    BiCG_Converged(row_no_prec) = ...
        R.BiCG_Converged;


    fast_Time(row_no_prec) = ...
        R.fast_Time;

    fast_RelRes(row_no_prec) = ...
        R.fast_RelRes;

    fast_OuterIter(row_no_prec) = ...
        R.fast_OuterIter;

    fast_Converged(row_no_prec) = ...
        R.fast_Converged;


    performance_time(row_no_prec,:) = ...
        R.performance_time;


    %% ====================================================================
    % CASE B:
    % ILU(0) RIGHT PRECONDITIONING
    %
    % We explicitly solve
    %
    %       A M^{-1} y = b
    %
    % where
    %
    %       M approximately L*U
    %
    % and recover
    %
    %       x = M^{-1} y.
    %
    % Therefore ALL THREE solvers operate on the SAME right-preconditioned
    % linear system.
    % =====================================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' CASE B: %s -- ILU(0) RIGHT PRECONDITIONING\n', ...
        matrix_names{p});
    fprintf('============================================================\n');


    fprintf('\nConstructing zero-fill ILU factors...\n');


    try

        tic;

        [L,U] = ilu(A,ilu_setup);

        ilu_elapsed = toc;


        ILU_Setup_Time(row_ilu) = ...
            ilu_elapsed;

        ILU_Status(row_ilu) = ...
            "OK";


        fprintf('ILU(0) setup completed.\n');

        fprintf('ILU setup time : %.4f s\n', ...
            ilu_elapsed);


        %% ---------------------------------------------------------------
        % Right-preconditioner application
        %
        %       M^{-1}v = U^{-1} L^{-1} v
        %
        % IMPORTANT:
        %
        % Do NOT explicitly form inv(L), inv(U), or inv(M).
        % ---------------------------------------------------------------

        Psolve = ...
            @(v) U\(L\v);


        %% ---------------------------------------------------------------
        % Right-preconditioned matrix operator:
        %
        %       A_prec y = A * M^{-1} y
        % ---------------------------------------------------------------

        Aprec = ...
            @(y) A*Psolve(y);


        %% ---------------------------------------------------------------
        % Since original x0 = 0, transformed y0 is also zero.
        % ---------------------------------------------------------------

        y0 = zeros(n,1);


        %% ---------------------------------------------------------------
        % Run ALL THREE solvers on the SAME transformed system
        % ---------------------------------------------------------------

        R = run_solver_case( ...
            Aprec, ...
            A, ...
            Psolve, ...
            b, ...
            y0, ...
            tol, ...
            restart, ...
            gmres_max_outer, ...
            bicg_max_iter, ...
            fast_max_outer, ...
            opts);


        GMRES_Time(row_ilu) = ...
            R.GMRES_Time;

        GMRES_RelRes(row_ilu) = ...
            R.GMRES_RelRes;

        GMRES_Steps(row_ilu) = ...
            R.GMRES_Steps;

        GMRES_Converged(row_ilu) = ...
            R.GMRES_Converged;


        BiCG_Time(row_ilu) = ...
            R.BiCG_Time;

        BiCG_RelRes(row_ilu) = ...
            R.BiCG_RelRes;

        BiCG_Iter(row_ilu) = ...
            R.BiCG_Iter;

        BiCG_Converged(row_ilu) = ...
            R.BiCG_Converged;


        fast_Time(row_ilu) = ...
            R.fast_Time;

        fast_RelRes(row_ilu) = ...
            R.fast_RelRes;

        fast_OuterIter(row_ilu) = ...
            R.fast_OuterIter;

        fast_Converged(row_ilu) = ...
            R.fast_Converged;


        performance_time(row_ilu,:) = ...
            R.performance_time;


        clear L U Psolve Aprec y0


    catch ME

        fprintf('\n');
        fprintf('ILU(0) PRECONDITIONING ERROR:\n');

        fprintf('%s\n', ...
            ME.message);

        fprintf('\n');

        fprintf( ...
            'The unpreconditioned case remains valid, but the\n');

        fprintf( ...
            'ILU-preconditioned case for this matrix will be marked INVALID.\n');


        ILU_Status(row_ilu) = ...
            "FAILED";

        ILU_Setup_Time(row_ilu) = ...
            NaN;


        % All iterative results remain failure values.

        GMRES_Time(row_ilu) = NaN;

        GMRES_RelRes(row_ilu) = Inf;

        GMRES_Converged(row_ilu) = false;


        BiCG_Time(row_ilu) = NaN;

        BiCG_RelRes(row_ilu) = Inf;

        BiCG_Converged(row_ilu) = false;


        fast_Time(row_ilu) = NaN;

        fast_RelRes(row_ilu) = Inf;

        fast_Converged(row_ilu) = false;


        performance_time(row_ilu,:) = Inf;

    end


    %% ====================================================================
    % Save progress immediately after each matrix
    % =====================================================================

    fprintf('\n');
    fprintf('Completed both cases for: %s\n', ...
        matrix_names{p});


    Results = table( ...
        Case, ...
        Matrix, ...
        Mode, ...
        N, ...
        NNZ, ...
        RHS_Source, ...
        ILU_Status, ...
        ILU_Setup_Time, ...
        GMRES_Time, ...
        GMRES_RelRes, ...
        GMRES_Steps, ...
        GMRES_Converged, ...
        BiCG_Time, ...
        BiCG_RelRes, ...
        BiCG_Iter, ...
        BiCG_Converged, ...
        fast_Time, ...
        fast_RelRes, ...
        fast_OuterIter, ...
        fast_Converged);


    writetable( ...
        Results, ...
        'PublicMatrix_20_CASE_ILU0_results.csv');


    save( ...
        'PublicMatrix_20_CASE_ILU0_results.mat', ...
        'Results', ...
        'performance_time', ...
        'tol', ...
        'restart', ...
        'gmres_max_outer', ...
        'bicg_max_iter', ...
        'fast_max_outer', ...
        'opts', ...
        'ilu_setup');


    %% --------------------------------------------------------------------
    % Release large matrix before loading the next problem
    % ---------------------------------------------------------------------

    clear A b x0 R

    drawnow;

    pause(1);

end


%% ========================================================================
% 12. Final results table
% ========================================================================

fprintf('\n\n');
fprintf('============================================================\n');
fprintf(' FINAL RESULTS: 20 BENCHMARK CASES\n');
fprintf('============================================================\n\n');


Results = table( ...
    Case, ...
    Matrix, ...
    Mode, ...
    N, ...
    NNZ, ...
    RHS_Source, ...
    ILU_Status, ...
    ILU_Setup_Time, ...
    GMRES_Time, ...
    GMRES_RelRes, ...
    GMRES_Steps, ...
    GMRES_Converged, ...
    BiCG_Time, ...
    BiCG_RelRes, ...
    BiCG_Iter, ...
    BiCG_Converged, ...
    fast_Time, ...
    fast_RelRes, ...
    fast_OuterIter, ...
    fast_Converged);


disp(Results);


%% ========================================================================
% 13. Save final data
% ========================================================================

writetable( ...
    Results, ...
    'PublicMatrix_20_CASE_ILU0_results.csv');


save( ...
    'PublicMatrix_20_CASE_ILU0_results.mat', ...
    'Results', ...
    'performance_time', ...
    'tol', ...
    'restart', ...
    'gmres_max_outer', ...
    'bicg_max_iter', ...
    'fast_max_outer', ...
    'opts', ...
    'ilu_setup');


%% ========================================================================
% 14. Row groups
% ========================================================================

rows_no_prec = ...
    1:2:numCases;

rows_ilu = ...
    2:2:numCases;


%% ========================================================================
% 15. Convergence summary: ALL 20 cases
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CONVERGENCE SUMMARY: ALL 20 CASES\n');
fprintf('============================================================\n\n');


gmres_success_all = ...
    sum(GMRES_Converged);

bicg_success_all = ...
    sum(BiCG_Converged);

fast_success_all = ...
    sum(fast_Converged);


fprintf( ...
    'GMRES(50) : %d / %d cases converged\n', ...
    gmres_success_all, ...
    numCases);

fprintf( ...
    'BiCGSTAB   : %d / %d cases converged\n', ...
    bicg_success_all, ...
    numCases);

fprintf( ...
    'fastGMRES  : %d / %d cases converged\n', ...
    fast_success_all, ...
    numCases);


%% ========================================================================
% 16. Convergence summary: WITHOUT preconditioning
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' WITHOUT PRECONDITIONING: 10 CASES\n');
fprintf('============================================================\n\n');


fprintf( ...
    'GMRES(50) : %d / 10\n', ...
    sum(GMRES_Converged(rows_no_prec)));

fprintf( ...
    'BiCGSTAB   : %d / 10\n', ...
    sum(BiCG_Converged(rows_no_prec)));

fprintf( ...
    'fastGMRES  : %d / 10\n', ...
    sum(fast_Converged(rows_no_prec)));


%% ========================================================================
% 17. Convergence summary: ILU(0)
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' ILU(0) RIGHT-PRECONDITIONED: 10 CASES\n');
fprintf('============================================================\n\n');


fprintf( ...
    'Successful ILU constructions : %d / 10\n\n', ...
    sum(ILU_Status(rows_ilu) == "OK"));


fprintf( ...
    'GMRES(50) : %d / 10\n', ...
    sum(GMRES_Converged(rows_ilu)));

fprintf( ...
    'BiCGSTAB   : %d / 10\n', ...
    sum(BiCG_Converged(rows_ilu)));

fprintf( ...
    'fastGMRES  : %d / 10\n', ...
    sum(fast_Converged(rows_ilu)));


%% ========================================================================
% 18. Performance profile input
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' PERFORMANCE PROFILE INPUT: 20 CASES\n');
fprintf('============================================================\n\n');

fprintf( ...
    'Inf = target tolerance %.1e was NOT achieved.\n\n', ...
    tol);


disp( ...
    array2table( ...
        performance_time, ...
        'VariableNames', ...
        {'GMRES50','BiCGSTAB','fastGMRES'}, ...
        'RowNames', ...
        cellstr(Case)));


%% ========================================================================
% 19. PERFORMANCE PROFILE:
% ALL 20 CASES
%
% PDF vector output.
% White background.
% ========================================================================

fig1 = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [200 150 900 600]);


pprofile_local( ...
    performance_time, ...
    methods);


title( ...
    'Runtime Performance Profile: 20 Benchmark Cases', ...
    'FontWeight', ...
    'normal');


exportgraphics( ...
    fig1, ...
    'Performance_Profile_20Cases.pdf', ...
    'ContentType', ...
    'vector');


%% ========================================================================
% 20. PERFORMANCE PROFILE:
% 10 UNPRECONDITIONED CASES
% ========================================================================

fig2 = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [200 150 900 600]);


pprofile_local( ...
    performance_time(rows_no_prec,:), ...
    methods);


title( ...
    'Runtime Performance Profile: 10 Unpreconditioned Cases', ...
    'FontWeight', ...
    'normal');


exportgraphics( ...
    fig2, ...
    'Performance_Profile_10_Unpreconditioned.pdf', ...
    'ContentType', ...
    'vector');


%% ========================================================================
% 21. PERFORMANCE PROFILE:
% 10 ILU(0)-PRECONDITIONED CASES
% ========================================================================

fig3 = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [200 150 900 600]);


pprofile_local( ...
    performance_time(rows_ilu,:), ...
    methods);


title( ...
    'Runtime Performance Profile: 10 ILU(0)-Preconditioned Cases', ...
    'FontWeight', ...
    'normal');


exportgraphics( ...
    fig3, ...
    'Performance_Profile_10_ILU0.pdf', ...
    'ContentType', ...
    'vector');


%% ========================================================================
% 22. Save plain-text summary
% ========================================================================

fid = fopen( ...
    'PublicMatrix_20_CASE_ILU0_summary.txt', ...
    'w');


fprintf(fid, ...
    'FINAL 20-CASE BENCHMARK\n\n');


fprintf(fid, ...
    '10 matrices tested without and with ILU(0) right preconditioning.\n\n');


fprintf(fid, ...
    'Target tolerance: %.1e\n\n', ...
    tol);


fprintf(fid, ...
    'ILU: zero-fill incomplete LU (ILU(0)).\n');

fprintf(fid, ...
    'ILU setup time is recorded separately and excluded from solver runtime.\n\n');


fprintf(fid, ...
    'ALL 20 CASES\n');

fprintf(fid, ...
    'GMRES(50): %d / %d converged\n', ...
    gmres_success_all, ...
    numCases);

fprintf(fid, ...
    'BiCGSTAB: %d / %d converged\n', ...
    bicg_success_all, ...
    numCases);

fprintf(fid, ...
    'fastGMRES: %d / %d converged\n\n', ...
    fast_success_all, ...
    numCases);


fprintf(fid, ...
    'WITHOUT PRECONDITIONING\n');

fprintf(fid, ...
    'GMRES(50): %d / 10\n', ...
    sum(GMRES_Converged(rows_no_prec)));

fprintf(fid, ...
    'BiCGSTAB: %d / 10\n', ...
    sum(BiCG_Converged(rows_no_prec)));

fprintf(fid, ...
    'fastGMRES: %d / 10\n\n', ...
    sum(fast_Converged(rows_no_prec)));


fprintf(fid, ...
    'WITH ILU(0) RIGHT PRECONDITIONING\n');

fprintf(fid, ...
    'ILU constructions successful: %d / 10\n', ...
    sum(ILU_Status(rows_ilu) == "OK"));

fprintf(fid, ...
    'GMRES(50): %d / 10\n', ...
    sum(GMRES_Converged(rows_ilu)));

fprintf(fid, ...
    'BiCGSTAB: %d / 10\n', ...
    sum(BiCG_Converged(rows_ilu)));

fprintf(fid, ...
    'fastGMRES: %d / 10\n', ...
    sum(fast_Converged(rows_ilu)));


fclose(fid);


%% ========================================================================
% 23. Final message
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' SAVED FILES\n');
fprintf('============================================================\n\n');


fprintf( ...
    'PublicMatrix_20_CASE_ILU0_results.csv\n');

fprintf( ...
    'PublicMatrix_20_CASE_ILU0_results.mat\n');

fprintf( ...
    'PublicMatrix_20_CASE_ILU0_summary.txt\n');

fprintf( ...
    'Performance_Profile_20Cases.pdf\n');

fprintf( ...
    'Performance_Profile_10_Unpreconditioned.pdf\n');

fprintf( ...
    'Performance_Profile_10_ILU0.pdf\n');


fprintf('\n');
fprintf('============================================================\n');
fprintf(' ALL 20 BENCHMARK CASES COMPLETED\n');
fprintf('============================================================\n');



%% ========================================================================
% LOCAL FUNCTION 1
%
% Run GMRES, BiCGSTAB and fastGMRES on ONE linear system.
%
%
% Aop:
%
%       Matrix or function handle used by the solver.
%
%
% Atrue:
%
%       Original physical matrix A.
%       Used to independently compute the TRUE residual.
%
%
% recover_x:
%
%       Maps solver variable back to original solution x.
%
%       Without preconditioning:
%
%               x = y
%
%       Right-preconditioning:
%
%               x = M^{-1} y
%
%
% Therefore convergence is ALWAYS judged using
%
%       ||b - A*x|| / ||b||
%
% for a fair comparison.
% ========================================================================

function R = run_solver_case( ...
    Aop, ...
    Atrue, ...
    recover_x, ...
    b, ...
    x0, ...
    tol, ...
    restart, ...
    gmres_max_outer, ...
    bicg_max_iter, ...
    fast_max_outer, ...
    opts)


    %% --------------------------------------------------------------------
    % Initialize result structure
    % ---------------------------------------------------------------------

    R.GMRES_Time = NaN;

    R.GMRES_RelRes = Inf;

    R.GMRES_Steps = NaN;

    R.GMRES_Converged = false;


    R.BiCG_Time = NaN;

    R.BiCG_RelRes = Inf;

    R.BiCG_Iter = NaN;

    R.BiCG_Converged = false;


    R.fast_Time = NaN;

    R.fast_RelRes = Inf;

    R.fast_OuterIter = NaN;

    R.fast_Converged = false;


    R.performance_time = ...
        Inf(1,3);


    %% ====================================================================
    % GMRES(50)
    % =====================================================================

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('GMRES(50)\n');
    fprintf('------------------------------------------------------------\n');


    try

        tic;


        [y,flag,~,iter,~] = ...
            gmres( ...
                Aop, ...
                b, ...
                restart, ...
                tol, ...
                gmres_max_outer, ...
                [], ...
                [], ...
                x0);


        % Convert transformed variable y to original solution x
        x = recover_x(y);


        elapsed = toc;


        %% ---------------------------------------------------------------
        % Independently evaluated TRUE original-system residual
        % ---------------------------------------------------------------

        true_relres = ...
            norm(b-Atrue*x)/norm(b);


        R.GMRES_Time = ...
            elapsed;

        R.GMRES_RelRes = ...
            true_relres;


        %% ---------------------------------------------------------------
        % MATLAB GMRES returns:
        %
        %       iter = [outer inner]
        %
        % Convert this into approximate total Krylov steps.
        % ---------------------------------------------------------------

        if numel(iter) == 2

            if iter(1) == 0

                steps = ...
                    iter(2);

            else

                steps = ...
                    (iter(1)-1)*restart + iter(2);

            end

            R.GMRES_Steps = ...
                steps;

        end


        R.GMRES_Converged = ...
            true_relres <= tol;


        if R.GMRES_Converged

            R.performance_time(1) = ...
                elapsed;

        end


        fprintf('MATLAB flag       : %d\n', ...
            flag);

        fprintf('Runtime           : %.4f s\n', ...
            elapsed);

        fprintf('Krylov steps      : %.0f\n', ...
            R.GMRES_Steps);

        fprintf('True relative res : %.8e\n', ...
            true_relres);


        if R.GMRES_Converged

            fprintf('STATUS            : CONVERGED\n');

        else

            fprintf('STATUS            : FAILED\n');

        end


        clear x y


    catch ME

        fprintf('\nGMRES ERROR:\n');

        fprintf('%s\n', ...
            ME.message);

    end


    %% ====================================================================
    % BiCGSTAB
    % =====================================================================

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('BiCGSTAB\n');
    fprintf('------------------------------------------------------------\n');


    try

        tic;


        [y,flag,~,iter,~] = ...
            bicgstab( ...
                Aop, ...
                b, ...
                tol, ...
                bicg_max_iter, ...
                [], ...
                [], ...
                x0);


        x = ...
            recover_x(y);


        elapsed = toc;


        true_relres = ...
            norm(b-Atrue*x)/norm(b);


        R.BiCG_Time = ...
            elapsed;

        R.BiCG_RelRes = ...
            true_relres;

        R.BiCG_Iter = ...
            iter;


        R.BiCG_Converged = ...
            true_relres <= tol;


        if R.BiCG_Converged

            R.performance_time(2) = ...
                elapsed;

        end


        fprintf('MATLAB flag       : %d\n', ...
            flag);

        fprintf('Runtime           : %.4f s\n', ...
            elapsed);

        fprintf('Iterations        : %.1f\n', ...
            iter);

        fprintf('True relative res : %.8e\n', ...
            true_relres);


        if R.BiCG_Converged

            fprintf('STATUS            : CONVERGED\n');

        else

            fprintf('STATUS            : FAILED\n');

        end


        clear x y


    catch ME

        fprintf('\nBiCGSTAB ERROR:\n');

        fprintf('%s\n', ...
            ME.message);

    end


    %% ====================================================================
    % fastGMRES
    % =====================================================================

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('fastGMRES\n');
    fprintf('------------------------------------------------------------\n');


    try

        %% ---------------------------------------------------------------
        % Identical MATLAB random state for every fastGMRES solve.
        %
        % This preserves reproducibility of the randomized sketch.
        % ---------------------------------------------------------------

        rng('default');


        tic;


        [y,resvec,restime] = ...
            fastgmres( ...
                Aop, ...
                b, ...
                tol, ...
                fast_max_outer, ...
                [], ...
                [], ...
                x0, ...
                opts);


        x = ...
            recover_x(y);


        elapsed = toc;


        true_relres = ...
            norm(b-Atrue*x)/norm(b);


        R.fast_Time = ...
            elapsed;

        R.fast_RelRes = ...
            true_relres;

        R.fast_OuterIter = ...
            length(resvec);


        R.fast_Converged = ...
            true_relres <= tol;


        if R.fast_Converged

            R.performance_time(3) = ...
                elapsed;

        end


        fprintf('Runtime           : %.4f s\n', ...
            elapsed);

        fprintf('Outer iterations  : %d\n', ...
            length(resvec));

        fprintf('True relative res : %.8e\n', ...
            true_relres);


        if R.fast_Converged

            fprintf('STATUS            : CONVERGED\n');

        else

            fprintf('STATUS            : FAILED\n');

        end


        clear x y resvec restime


    catch ME

        fprintf('\nfastGMRES ERROR:\n');

        fprintf('%s\n', ...
            ME.message);


        if contains( ...
                lower(ME.message), ...
                'memory')

            fprintf('\nOUT OF MEMORY detected.\n');

        end

    end

end



%% ========================================================================
% LOCAL FUNCTION 2
%
% Load SuiteSparse / ordinary MAT file.
%
% IMPORTANT:
%
% Supplied right-hand sides are intentionally ignored.
%
% Every matrix uses:
%
%       x_true = reproducible random vector
%       b      = A*x_true
%
% This preserves the SAME experimental protocol used in the previous
% 10-matrix benchmark.
% ========================================================================

function [A,b,rhs_source] = ...
    load_matrix_problem(filename,rhs_seed)


    S = load(filename);


    A = [];


    %% --------------------------------------------------------------------
    % SuiteSparse format:
    %
    %       Problem.A
    % ---------------------------------------------------------------------

    if isfield(S,'Problem')

        P = S.Problem;


        if isfield(P,'A')

            A = P.A;

        end

    end


    %% --------------------------------------------------------------------
    % Ordinary top-level format:
    %
    %       A
    % ---------------------------------------------------------------------

    if isempty(A) && isfield(S,'A')

        A = S.A;

    end


    %% --------------------------------------------------------------------
    % Automatically search for largest square numeric matrix
    % ---------------------------------------------------------------------

    if isempty(A)

        names = fieldnames(S);

        best_n = 0;


        for j = 1:length(names)

            X = S.(names{j});


            if isnumeric(X) && ...
                    ismatrix(X)

                [m,n] = size(X);


                if m == n && ...
                        m > best_n && ...
                        m > 1

                    A = X;

                    best_n = m;

                end

            end

        end

    end


    %% --------------------------------------------------------------------
    % Matrix must have been found
    % ---------------------------------------------------------------------

    if isempty(A)

        error( ...
            'No square matrix was found inside %s.', ...
            filename);

    end


    %% --------------------------------------------------------------------
    % Convert to sparse
    % ---------------------------------------------------------------------

    A = sparse(A);


    %% --------------------------------------------------------------------
    % Reproducible manufactured RHS
    % ---------------------------------------------------------------------

    rng(rhs_seed,'twister');


    x_true = ...
        randn(size(A,1),1);


    b = ...
        A*x_true;


    b = ...
        full(b(:));


    clear x_true


    rhs_source = ...
        "Manufactured: b = A*x_true";

end



%% ========================================================================
% LOCAL FUNCTION 3
%
% Runtime Performance Profile
%
% Smaller runtime is better.
%
% For method j:
%
%       P_j(alpha)
%
% is the fraction of benchmark problems where method j takes no more than
% alpha times the fastest successful solver.
%
%
% IMPORTANT:
%
% If NO solver reaches tolerance for a given problem, there is no finite
% reference runtime for that problem. Such rows are excluded from the
% performance-ratio denominator.
%
% A solver failure on a problem solved by another method remains penalized.
% ========================================================================

function pprofile_local(M,methods)


    %% --------------------------------------------------------------------
    % Identify benchmark rows solved by at least one method
    % ---------------------------------------------------------------------

    valid_problem = ...
        any(isfinite(M),2);


    Mvalid = ...
        M(valid_problem,:);


    nProblems = ...
        size(Mvalid,1);


    nMethods = ...
        size(Mvalid,2);


    if nProblems == 0

        warning( ...
            'No benchmark problem was solved by any method.');

        text( ...
            0.5, ...
            0.5, ...
            'No successful benchmark cases', ...
            'HorizontalAlignment', ...
            'center');

        axis off;

        return;

    end


    excluded = ...
        size(M,1) - nProblems;


    fprintf('\nPerformance profile:\n');

    fprintf( ...
        '  Included cases = %d\n', ...
        nProblems);

    fprintf( ...
        '  Excluded cases = %d (no method converged)\n', ...
        excluded);


    %% --------------------------------------------------------------------
    % Alpha range
    % ---------------------------------------------------------------------

    alpha = ...
        linspace(1,4,200);


    P = ...
        zeros(length(alpha),nMethods);


    %% --------------------------------------------------------------------
    % Fastest successful runtime for each problem
    % ---------------------------------------------------------------------

    best = ...
        min(Mvalid,[],2);


    %% --------------------------------------------------------------------
    % Performance profile
    % ---------------------------------------------------------------------

    for k = 1:length(alpha)

        for j = 1:nMethods


            successful = ...
                isfinite(Mvalid(:,j)) & ...
                Mvalid(:,j) <= ...
                alpha(k).*best;


            P(k,j) = ...
                sum(successful)/nProblems;

        end

    end


    %% --------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------

    hold on;


    marker = {
        '-o'
        '-s'
        '-^'
    };


    for j = 1:nMethods

        plot( ...
            alpha, ...
            P(:,j), ...
            marker{j}, ...
            'LineWidth', ...
            1.6, ...
            'MarkerIndices', ...
            1:20:length(alpha));

    end


    grid on;

    box on;


    axis([1 4 0 1.05]);


    xlabel( ...
        'Performance ratio \alpha');


    ylabel( ...
        'Fraction of test problems');


    %% --------------------------------------------------------------------
    % Failure counts
    %
    % Only count failures on VALID benchmark rows.
    % ---------------------------------------------------------------------

    failures = ...
        sum(~isfinite(Mvalid),1);


    leg = ...
        cell(1,nMethods);


    for j = 1:nMethods

        if failures(j) == 0

            leg{j} = ...
                methods{j};

        else

            leg{j} = ...
                sprintf( ...
                    '%s (%d failures)', ...
                    methods{j}, ...
                    failures(j));

        end

    end


    legend( ...
        leg, ...
        'Location', ...
        'southeast');

end