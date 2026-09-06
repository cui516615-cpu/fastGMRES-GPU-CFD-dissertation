%% ========================================================================
% Run_PublicMatrix_10_FINAL.m
%
% Final 10-matrix local-PC benchmark for Chapter 5
%
% Matrices:
%   1.  naca12A.mat
%   2.  ML_Laplace.mat
%   3.  vas_stokes_1M.mat
%   4.  Transport.mat
%   5.  cavity10.mat
%   6.  ns3Da.mat
%   7.  atmosmodl.mat
%   8.  atmosmodd.mat
%   9.  GT01R.mat
%   10. poisson3Da.mat
%
% Solvers:
%   GMRES(50)
%   BiCGSTAB
%   fastGMRES
%
% Common target tolerance:
%   1e-7
%
% All public matrices use a reproducible manufactured RHS:
%
%       b = A*x_true
%
% Designed for a 16 GB RAM laptop.
% ========================================================================

clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf(' FINAL 10-MATRIX BENCHMARK FOR CHAPTER 5\n');
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

numProblems = length(files);


%% ========================================================================
% 2. Check that all matrix files exist
% ========================================================================

fprintf('Checking matrix files...\n');

for p = 1:numProblems

    if ~isfile(files{p})

        error( ...
            'Matrix file "%s" was not found in the current MATLAB folder.', ...
            files{p});

    end

    fprintf('  [%2d/10] %-22s : OK\n', p, files{p});

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

% 40 restart cycles:
%
%       40 * 50 = 2000 maximum Krylov steps
%
gmres_max_outer = 40;


% -------------------------------------------------------------------------
% BiCGSTAB
% -------------------------------------------------------------------------

% Approximately two SpMVs per iteration
bicg_max_iter = 1000;


% -------------------------------------------------------------------------
% fastGMRES
% -------------------------------------------------------------------------
%
% These parameters are kept identical for all ten public matrices.
%
% m = 80 and maxit = 150 are chosen so that the million-scale matrices
% remain feasible on the 16 GB test computer.
%

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

fprintf('Target tolerance        = %.1e\n', tol);

fprintf('\nGMRES restart           = %d\n', restart);
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
% 4. Check fastgmres.m
% ========================================================================

if exist('fastgmres', 'file') ~= 2

    error( ...
        'fastgmres.m was not found in the current MATLAB folder or path.');

end


%% ========================================================================
% 5. Allocate result arrays
% ========================================================================

N = zeros(numProblems,1);

NNZ = zeros(numProblems,1);

RHS_Source = strings(numProblems,1);


% -------------------------------------------------------------------------
% GMRES results
% -------------------------------------------------------------------------

GMRES_Time = NaN(numProblems,1);

GMRES_RelRes = Inf(numProblems,1);

GMRES_Steps = NaN(numProblems,1);

GMRES_Converged = false(numProblems,1);


% -------------------------------------------------------------------------
% BiCGSTAB results
% -------------------------------------------------------------------------

BiCG_Time = NaN(numProblems,1);

BiCG_RelRes = Inf(numProblems,1);

BiCG_Iter = NaN(numProblems,1);

BiCG_Converged = false(numProblems,1);


% -------------------------------------------------------------------------
% fastGMRES results
% -------------------------------------------------------------------------

fast_Time = NaN(numProblems,1);

fast_RelRes = Inf(numProblems,1);

fast_OuterIter = NaN(numProblems,1);

fast_Converged = false(numProblems,1);


% -------------------------------------------------------------------------
% Runtime matrix used by performance profile
%
% Column 1 = GMRES(50)
% Column 2 = BiCGSTAB
% Column 3 = fastGMRES
%
% Inf = target tolerance was not achieved
% -------------------------------------------------------------------------

performance_time = Inf(numProblems,3);


%% ========================================================================
% 6. Main benchmark loop
% ========================================================================

for p = 1:numProblems

    fprintf('\n\n');
    fprintf('============================================================\n');
    fprintf(' MATRIX %d / %d : %s\n', ...
        p, numProblems, matrix_names{p});
    fprintf('============================================================\n');


    %% --------------------------------------------------------------------
    % Load matrix
    % ---------------------------------------------------------------------

    fprintf('\n[1] Loading matrix...\n');

    [A,b,rhs_source] = ...
        load_matrix_problem(files{p},1000+p);

    A = sparse(A);

    b = full(b(:));

    n = size(A,1);


    % ---------------------------------------------------------------------
    % Basic checks
    % ---------------------------------------------------------------------

    if size(A,1) ~= size(A,2)

        error('%s is not square.',files{p});

    end

    if length(b) ~= n

        error('RHS dimension mismatch for %s.',files{p});

    end

    if norm(b) == 0

        error('RHS is zero for %s.',files{p});

    end


    % ---------------------------------------------------------------------
    % Store matrix information
    % ---------------------------------------------------------------------

    N(p) = n;

    NNZ(p) = nnz(A);

    RHS_Source(p) = rhs_source;


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


    % Common initial guess
    x0 = zeros(n,1);


    %% ====================================================================
    % 2. GMRES(50)
    % =====================================================================

    fprintf('\n------------------------------------------------------------\n');
    fprintf('[2] GMRES(50)\n');
    fprintf('------------------------------------------------------------\n');

    try

        tic;

        [x,flag,~,iter,~] = ...
            gmres( ...
                A, ...
                b, ...
                restart, ...
                tol, ...
                gmres_max_outer, ...
                [], ...
                [], ...
                x0);

        elapsed = toc;


        % -----------------------------------------------------------------
        % Independently evaluated TRUE relative residual
        % -----------------------------------------------------------------

        true_relres = ...
            norm(b-A*x)/norm(b);


        GMRES_Time(p) = elapsed;

        GMRES_RelRes(p) = true_relres;


        % -----------------------------------------------------------------
        % Convert MATLAB [outer inner] iteration output
        % into approximate total Krylov steps
        % -----------------------------------------------------------------

        if numel(iter) == 2

            if iter(1) == 0

                steps = iter(2);

            else

                steps = ...
                    (iter(1)-1)*restart + iter(2);

            end

            GMRES_Steps(p) = steps;

        end


        GMRES_Converged(p) = ...
            true_relres <= tol;


        if GMRES_Converged(p)

            performance_time(p,1) = elapsed;

        end


        fprintf('MATLAB flag       : %d\n',flag);

        fprintf('Runtime           : %.4f s\n',elapsed);

        fprintf('Krylov steps      : %.0f\n', ...
            GMRES_Steps(p));

        fprintf('True relative res : %.8e\n', ...
            true_relres);


        if GMRES_Converged(p)

            fprintf('STATUS            : CONVERGED\n');

        else

            fprintf('STATUS            : FAILED\n');

        end


        clear x


    catch ME

        fprintf('\nGMRES ERROR:\n');

        fprintf('%s\n',ME.message);

        GMRES_Time(p) = NaN;

        GMRES_RelRes(p) = Inf;

        GMRES_Converged(p) = false;

    end


    %% ====================================================================
    % 3. BiCGSTAB
    % =====================================================================

    fprintf('\n------------------------------------------------------------\n');
    fprintf('[3] BiCGSTAB\n');
    fprintf('------------------------------------------------------------\n');

    try

        tic;

        [x,flag,~,iter,~] = ...
            bicgstab( ...
                A, ...
                b, ...
                tol, ...
                bicg_max_iter, ...
                [], ...
                [], ...
                x0);

        elapsed = toc;


        % -----------------------------------------------------------------
        % True relative residual
        % -----------------------------------------------------------------

        true_relres = ...
            norm(b-A*x)/norm(b);


        BiCG_Time(p) = elapsed;

        BiCG_RelRes(p) = true_relres;

        BiCG_Iter(p) = iter;


        BiCG_Converged(p) = ...
            true_relres <= tol;


        if BiCG_Converged(p)

            performance_time(p,2) = elapsed;

        end


        fprintf('MATLAB flag       : %d\n',flag);

        fprintf('Runtime           : %.4f s\n',elapsed);

        fprintf('Iterations        : %.1f\n',iter);

        fprintf('True relative res : %.8e\n', ...
            true_relres);


        if BiCG_Converged(p)

            fprintf('STATUS            : CONVERGED\n');

        else

            fprintf('STATUS            : FAILED\n');

        end


        clear x


    catch ME

        fprintf('\nBiCGSTAB ERROR:\n');

        fprintf('%s\n',ME.message);

        BiCG_Time(p) = NaN;

        BiCG_RelRes(p) = Inf;

        BiCG_Converged(p) = false;

    end


    %% ====================================================================
    % 4. fastGMRES
    % =====================================================================

    fprintf('\n------------------------------------------------------------\n');
    fprintf('[4] fastGMRES\n');
    fprintf('------------------------------------------------------------\n');

    try

        % ---------------------------------------------------------------
        % Same random sketch state for every test problem.
        %
        % This makes the benchmark reproducible.
        % ---------------------------------------------------------------

        rng('default');


        tic;

        [x,resvec,restime] = ...
            fastgmres( ...
                A, ...
                b, ...
                tol, ...
                fast_max_outer, ...
                [], ...
                [], ...
                x0, ...
                opts);

        elapsed = toc;


        % -----------------------------------------------------------------
        % Independently evaluated TRUE relative residual
        % -----------------------------------------------------------------

        true_relres = ...
            norm(b-A*x)/norm(b);


        fast_Time(p) = elapsed;

        fast_RelRes(p) = true_relres;

        fast_OuterIter(p) = length(resvec);


        fast_Converged(p) = ...
            true_relres <= tol;


        if fast_Converged(p)

            performance_time(p,3) = elapsed;

        end


        fprintf('Runtime           : %.4f s\n',elapsed);

        fprintf('Outer iterations  : %d\n', ...
            length(resvec));

        fprintf('True relative res : %.8e\n', ...
            true_relres);


        if fast_Converged(p)

            fprintf('STATUS            : CONVERGED\n');

        else

            fprintf('STATUS            : FAILED\n');

        end


        clear x resvec restime


    catch ME

        fprintf('\nfastGMRES ERROR:\n');

        fprintf('%s\n',ME.message);


        if contains(lower(ME.message),'memory')

            fprintf('\nOUT OF MEMORY detected.\n');

        end


        fast_Time(p) = NaN;

        fast_RelRes(p) = Inf;

        fast_Converged(p) = false;

    end


    %% ====================================================================
    % 5. Save progress immediately after every matrix
    % =====================================================================

    fprintf('\nCompleted: %s\n', ...
        matrix_names{p});


    Matrix = string(matrix_names);


    Results = table( ...
        Matrix, ...
        N, ...
        NNZ, ...
        RHS_Source, ...
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
        'PublicMatrix_10_FINAL_results.csv');


    save( ...
        'PublicMatrix_10_FINAL_results.mat', ...
        'Results', ...
        'performance_time', ...
        'tol', ...
        'restart', ...
        'gmres_max_outer', ...
        'bicg_max_iter', ...
        'fast_max_outer', ...
        'opts');


    % ---------------------------------------------------------------------
    % Release large matrix before loading next problem
    % ---------------------------------------------------------------------

    clear A b x0

    drawnow;

    pause(1);

end


%% ========================================================================
% 7. Final results table
% ========================================================================

fprintf('\n\n============================================================\n');
fprintf(' FINAL RESULTS: 10 MATRICES\n');
fprintf('============================================================\n\n');


Matrix = string(matrix_names);


Results = table( ...
    Matrix, ...
    N, ...
    NNZ, ...
    RHS_Source, ...
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
% 8. Save final results
% ========================================================================

writetable( ...
    Results, ...
    'PublicMatrix_10_FINAL_results.csv');


save( ...
    'PublicMatrix_10_FINAL_results.mat', ...
    'Results', ...
    'performance_time', ...
    'tol', ...
    'restart', ...
    'gmres_max_outer', ...
    'bicg_max_iter', ...
    'fast_max_outer', ...
    'opts');


%% ========================================================================
% 9. Print convergence summary
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' CONVERGENCE SUMMARY\n');
fprintf('============================================================\n\n');


gmres_success = ...
    sum(GMRES_Converged);

bicg_success = ...
    sum(BiCG_Converged);

fast_success = ...
    sum(fast_Converged);


fprintf( ...
    'GMRES(50) : %d / %d matrices converged\n', ...
    gmres_success, ...
    numProblems);

fprintf( ...
    'BiCGSTAB   : %d / %d matrices converged\n', ...
    bicg_success, ...
    numProblems);

fprintf( ...
    'fastGMRES  : %d / %d matrices converged\n', ...
    fast_success, ...
    numProblems);


fprintf('\nSuccess rates:\n');

fprintf( ...
    'GMRES(50) : %.1f %%\n', ...
    100*gmres_success/numProblems);

fprintf( ...
    'BiCGSTAB   : %.1f %%\n', ...
    100*bicg_success/numProblems);

fprintf( ...
    'fastGMRES  : %.1f %%\n', ...
    100*fast_success/numProblems);


%% ========================================================================
% 10. Performance profile input
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' PERFORMANCE PROFILE INPUT\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('Inf = target tolerance 1e-7 was NOT achieved.\n\n');


disp( ...
    array2table( ...
        performance_time, ...
        'VariableNames', ...
        {'GMRES50','BiCGSTAB','fastGMRES'}, ...
        'RowNames', ...
        matrix_names));


%% ========================================================================
% 11. Plot performance profile
% ========================================================================

figure( ...
    'Position', ...
    [200 150 900 600]);


pprofile_local( ...
    performance_time, ...
    methods);


title( ...
    'Runtime Performance Profile: 10 Non-Symmetric Matrices');


exportgraphics( ...
    gcf, ...
    'Performance_Profile_10Matrices.png', ...
    'Resolution', ...
    300);


%% ========================================================================
% 12. Save plain-text summary
% ========================================================================

fid = fopen( ...
    'PublicMatrix_10_FINAL_summary.txt', ...
    'w');


fprintf(fid, ...
    'FINAL 10-MATRIX BENCHMARK\n\n');

fprintf(fid, ...
    'Tolerance: %.1e\n\n', ...
    tol);


fprintf(fid, ...
    'GMRES(50): %d / %d converged (%.1f%%)\n', ...
    gmres_success, ...
    numProblems, ...
    100*gmres_success/numProblems);


fprintf(fid, ...
    'BiCGSTAB: %d / %d converged (%.1f%%)\n', ...
    bicg_success, ...
    numProblems, ...
    100*bicg_success/numProblems);


fprintf(fid, ...
    'fastGMRES: %d / %d converged (%.1f%%)\n', ...
    fast_success, ...
    numProblems, ...
    100*fast_success/numProblems);


fclose(fid);


%% ========================================================================
% 13. Final message
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' SAVED FILES\n');
fprintf('============================================================\n\n');


fprintf( ...
    'PublicMatrix_10_FINAL_results.csv\n');

fprintf( ...
    'PublicMatrix_10_FINAL_results.mat\n');

fprintf( ...
    'PublicMatrix_10_FINAL_summary.txt\n');

fprintf( ...
    'Performance_Profile_10Matrices.png\n');


fprintf('\n============================================================\n');
fprintf(' ALL 10 TESTS COMPLETED\n');
fprintf('============================================================\n');



%% ========================================================================
% LOCAL FUNCTION 1
%
% Load SuiteSparse / ordinary MAT file.
%
% IMPORTANT:
%
% For experimental consistency, supplied right-hand sides are intentionally
% ignored in the 10-matrix public benchmark.
%
% Every matrix uses:
%
%       x_true = reproducible random vector
%       b      = A*x_true
%
% Therefore all ten public test problems use the same RHS construction rule.
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
    % Automatically search for the largest square numeric matrix.
    %
    % This also handles files such as naca12A.mat.
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
    % Convert to sparse representation
    % ---------------------------------------------------------------------

    A = sparse(A);


    %% --------------------------------------------------------------------
    % Generate reproducible manufactured RHS
    %
    % Each matrix receives a different, but fixed, random seed.
    % ---------------------------------------------------------------------

    rng(rhs_seed,'twister');


    x_true = ...
        randn(size(A,1),1);


    b = ...
        A*x_true;


    b = full(b(:));


    clear x_true


    rhs_source = ...
        "Manufactured: b = A*x_true";

end



%% ========================================================================
% LOCAL FUNCTION 2
%
% Runtime Performance Profile
%
% Smaller runtime is better.
%
% A solver that does not reach the target residual tolerance receives:
%
%       Inf
%
% For alpha >= 1:
%
% P_j(alpha)
%
% is the fraction of benchmark problems on which method j has runtime
% no more than alpha times the fastest successful method.
% ========================================================================

function pprofile_local(M,methods)


    alpha = ...
        linspace(1,4,200);


    nProblems = ...
        size(M,1);


    nMethods = ...
        size(M,2);


    P = ...
        zeros(length(alpha),nMethods);


    %% --------------------------------------------------------------------
    % Best successful runtime for every matrix
    % ---------------------------------------------------------------------

    best = ...
        min(M,[],2);


    %% --------------------------------------------------------------------
    % Construct profile
    % ---------------------------------------------------------------------

    for k = 1:length(alpha)

        for j = 1:nMethods


            successful = ...
                isfinite(M(:,j)) & ...
                isfinite(best) & ...
                M(:,j) <= ...
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
            1.5, ...
            'MarkerIndices', ...
            1:20:length(alpha));

    end


    grid on;

    box on;


    axis([1 4 0 1.05]);


    xlabel('\alpha');

    ylabel('Fraction of test problems');


    %% --------------------------------------------------------------------
    % Failure counts
    % ---------------------------------------------------------------------

    failures = ...
        sum(~isfinite(M),1);


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