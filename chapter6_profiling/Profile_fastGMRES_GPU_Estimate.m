% =========================================================================
% Script: Profile_fastGMRES_GPU_Estimate.m
%
% Purpose:
%   Profile the SAME FluAtm fastGMRES case used in Chapter 5 in order to
%   construct an evidence-based estimate of future GPU performance.
%
% IMPORTANT:
%
%   This script does NOT change the numerical experiment.
%
%   Matrix:
%       FluAtm momentum matrix
%
%   Solver:
%       fastGMRES
%
%   Tolerance:
%       1e-7
%
%   Solver call:
%       fastgmres(A,b,tol)
%
% The script performs:
%
%   RUN 1:
%       Normal fastGMRES execution.
%       Used for the TRUE baseline CPU runtime.
%
%   MICROBENCHMARKS:
%       Measure representative CPU costs of
%           - sparse matrix-vector multiplication (SpMV)
%           - dot product
%           - AXPY
%           - vector norm
%           - vector scaling
%
%   RUN 2:
%       fastGMRES with MATLAB Profiler enabled.
%       Used ONLY to identify expensive functions and code lines.
%
% NO figures or files are generated.
%
% =========================================================================

clear;
clc;
close all;


fprintf('\n');
fprintf('============================================================\n');
fprintf(' FASTGMRES GPU-PERFORMANCE PROFILING EXPERIMENT\n');
fprintf(' FluAtm Momentum Matrix\n');
fprintf('============================================================\n\n');


%% ========================================================================
% 1. Check required files
% ========================================================================

row_file = 'fluatm_matrix_row.dat';
col_file = 'fluatm_matrix_col.dat';
val_file = 'fluatm_matrix_val.dat';
rhs_file = 'fluatm_rhs.dat';


files = {
    row_file
    col_file
    val_file
    rhs_file
};


fprintf('[Phase 1] Checking input files...\n');


for k = 1:length(files)

    if ~isfile(files{k})

        error( ...
            'Required input file "%s" was not found.', ...
            files{k});

    end

    fprintf('  %-28s : OK\n',files{k});

end


if exist('fastgmres','file') ~= 2

    error( ...
        ['fastgmres.m was not found in the current folder ' ...
         'or MATLAB path.']);

end


fprintf('  %-28s : OK\n\n','fastgmres.m');


%% ========================================================================
% 2. Read CSR row pointers
% ========================================================================

fprintf('[Phase 1] Reading FluAtm matrix...\n');


fid = fopen(row_file,'rb');

if fid == -1
    error('Failed to open %s.',row_file);
end


n_row_ptr = fread(fid,1,'int32');

row_ptr = fread(fid,n_row_ptr,'int32');

fclose(fid);


if length(row_ptr) ~= n_row_ptr

    error('Unexpected end of row-pointer file.');

end


nrows = n_row_ptr - 1;


%% ========================================================================
% 3. Read column indices
% ========================================================================

fid = fopen(col_file,'rb');

if fid == -1
    error('Failed to open %s.',col_file);
end


nnz_col = fread(fid,1,'int32');

col_idx = fread(fid,nnz_col,'int32');

fclose(fid);


if length(col_idx) ~= nnz_col

    error('Unexpected end of column-index file.');

end


%% ========================================================================
% 4. Read numerical values
% ========================================================================

fid = fopen(val_file,'rb');

if fid == -1
    error('Failed to open %s.',val_file);
end


nnz_val = fread(fid,1,'int32');

values = fread(fid,nnz_val,'double');

fclose(fid);


if length(values) ~= nnz_val

    error('Unexpected end of matrix-value file.');

end


%% ========================================================================
% 5. Read RHS
% ========================================================================

fid = fopen(rhs_file,'rb');

if fid == -1
    error('Failed to open %s.',rhs_file);
end


rhs_len = fread(fid,1,'int32');

b = fread(fid,rhs_len,'double');

fclose(fid);


b = full(b(:));


%% ========================================================================
% 6. Consistency checks
% ========================================================================

if nnz_col ~= nnz_val

    error( ...
        ['Number of column indices (%d) differs from ' ...
         'number of matrix values (%d).'], ...
        nnz_col, ...
        nnz_val);

end


if rhs_len ~= nrows

    error( ...
        'RHS dimension (%d) differs from matrix dimension (%d).', ...
        rhs_len, ...
        nrows);

end


%% ========================================================================
% 7. Detect column-index convention
% ========================================================================

if min(col_idx) == 0

    fprintf('  Detected 0-based column indices.\n');

    col_idx = col_idx + 1;

elseif min(col_idx) >= 1

    fprintf('  Column indices already appear to be 1-based.\n');

else

    error('Negative column index detected.');

end


%% ========================================================================
% 8. Convert CSR row pointers
% ========================================================================

row_counts = diff(row_ptr);


if any(row_counts < 0)

    error('CSR row pointers are not monotonically increasing.');

end


if sum(double(row_counts)) ~= nnz_val

    error( ...
        ['CSR row pointers predict %d nonzeros, ' ...
         'whereas the value file contains %d.'], ...
        sum(double(row_counts)), ...
        nnz_val);

end


row_idx = ...
    repelem( ...
        (1:nrows)', ...
        double(row_counts));


%% ========================================================================
% 9. Assemble sparse matrix
% ========================================================================

A = sparse( ...
    double(row_idx), ...
    double(col_idx), ...
    values, ...
    nrows, ...
    nrows);


norm_b = norm(b);


if norm_b == 0

    error('RHS vector has zero norm.');

end


fprintf('\n');
fprintf('============================================================\n');
fprintf(' MATRIX INFORMATION\n');
fprintf('============================================================\n');

fprintf('Dimension        : %d x %d\n', ...
    nrows,nrows);

fprintf('NNZ              : %d\n', ...
    nnz(A));

fprintf('Density          : %.6e\n', ...
    nnz(A)/(double(nrows)*double(nrows)));

fprintf('RHS norm         : %.8e\n', ...
    norm_b);

fprintf('============================================================\n\n');


% Free unnecessary CSR arrays
clear row_ptr row_counts row_idx col_idx values


%% ========================================================================
% 10. Solver configuration
%
% Keep EXACTLY the same tolerance used in the Chapter 5 FluAtm case.
% ========================================================================

tol = 1e-7;


fprintf('Solver           : fastGMRES\n');

fprintf('Tolerance        : %.1e\n',tol);

fprintf('Matrix           : FluAtm momentum matrix\n\n');


%% ========================================================================
% 11. Warm up basic MATLAB sparse operation
%
% This is NOT a solver warm-up.
%
% It only avoids including first-use overhead in the later kernel
% microbenchmark.
% ========================================================================

fprintf('[Phase 2] Warming up basic sparse operations...\n');


rng('default');


v_warm = randn(nrows,1);

tmp_warm = A*v_warm;


clear tmp_warm v_warm


fprintf('  Warm-up complete.\n\n');


%% ========================================================================
% 12. RUN 1:
%     NORMAL fastGMRES BASELINE
%
% This is the important runtime measurement.
%
% The MATLAB profiler is OFF here.
% ========================================================================

fprintf('============================================================\n');
fprintf(' RUN 1: NORMAL FASTGMRES BASELINE\n');
fprintf('============================================================\n\n');


profile clear;
profile off;


global TOTAL_SPMV_COUNT;

TOTAL_SPMV_COUNT = 0;


% Use reproducible randomized sketch
rng('default');


warning_state = warning;

warning('off','all');


fprintf('Running fastGMRES without profiler...\n');


tic;


[x_fast, ...
 rv_fast, ...
 rt_fast, ...
 spmv_history] = ...
    fastgmres( ...
        A, ...
        b, ...
        tol);


baseline_time = toc;


warning(warning_state);


%% ------------------------------------------------------------------------
% True residual
% -------------------------------------------------------------------------

true_relres = ...
    norm(b-A*x_fast)/norm_b;


%% ------------------------------------------------------------------------
% Recover SpMV information
% -------------------------------------------------------------------------

total_spmv = ...
    double(TOTAL_SPMV_COUNT);


% If the global counter is unavailable for any reason,
% recover the final cumulative count from spmv_history.
if total_spmv <= 0 && ~isempty(spmv_history)

    total_spmv = ...
        double(spmv_history(end));

end


% The previous Chapter 5 convention defines:
%
% number of outer iterations = length(rv_fast)
%
outer_iterations = ...
    length(rv_fast);


% One outer matrix-vector application per outer FGMRES step.
%
% The remainder therefore estimates cumulative inner-sGMRES SpMVs.
inner_spmv = ...
    max( ...
        total_spmv - outer_iterations, ...
        0);


if outer_iterations > 0

    average_inner_steps = ...
        inner_spmv / outer_iterations;

else

    average_inner_steps = NaN;

end


fprintf('\n');
fprintf('---------------- BASELINE RESULT ----------------\n');

fprintf('CPU fastGMRES runtime       : %.6f s\n', ...
    baseline_time);

fprintf('True relative residual      : %.8e\n', ...
    true_relres);

fprintf('Outer iterations            : %d\n', ...
    outer_iterations);

fprintf('Total recorded SpMVs        : %.0f\n', ...
    total_spmv);

fprintf('Estimated inner SpMVs       : %.0f\n', ...
    inner_spmv);

fprintf('Average inner steps / outer : %.3f\n', ...
    average_inner_steps);

fprintf('-------------------------------------------------\n\n');


%% ========================================================================
% 13. CPU KERNEL MICROBENCHMARKS
%
% We now measure representative costs of operations used by Krylov methods.
%
% These timings DO NOT by themselves represent the exact fastGMRES
% component times. They are used later to build an approximate performance
% model.
% ========================================================================

fprintf('============================================================\n');
fprintf(' CPU KERNEL MICROBENCHMARKS\n');
fprintf('============================================================\n\n');


rng(12345,'twister');


v1 = randn(nrows,1);

v2 = randn(nrows,1);

alpha_scalar = 0.713;


%% ------------------------------------------------------------------------
% SpMV benchmark
%
% 8 blocks x 10 operations
% -------------------------------------------------------------------------

fprintf('Benchmarking sparse matrix-vector multiplication...\n');


spmv_time = ...
    benchmark_spmv( ...
        A, ...
        v1, ...
        8, ...
        10);


%% ------------------------------------------------------------------------
% Dot product
%
% 8 blocks x 100 operations
% -------------------------------------------------------------------------

fprintf('Benchmarking dot products...\n');


dot_time = ...
    benchmark_dot( ...
        v1, ...
        v2, ...
        8, ...
        100);


%% ------------------------------------------------------------------------
% AXPY
%
% z = x + alpha*y
% -------------------------------------------------------------------------

fprintf('Benchmarking AXPY-style vector updates...\n');


axpy_time = ...
    benchmark_axpy( ...
        v1, ...
        v2, ...
        alpha_scalar, ...
        8, ...
        100);


%% ------------------------------------------------------------------------
% Norm
% -------------------------------------------------------------------------

fprintf('Benchmarking vector norms...\n');


norm_time = ...
    benchmark_norm( ...
        v1, ...
        8, ...
        100);


%% ------------------------------------------------------------------------
% Scaling
% -------------------------------------------------------------------------

fprintf('Benchmarking vector scaling...\n');


scale_time = ...
    benchmark_scale( ...
        v1, ...
        alpha_scalar, ...
        8, ...
        100);


%% ========================================================================
% 14. Estimate cumulative SpMV contribution
%
% IMPORTANT:
%
% This is an independent microbenchmark-based estimate.
% It is NOT a direct internal fastGMRES timer.
% ========================================================================

estimated_spmv_total = ...
    total_spmv * spmv_time;


estimated_spmv_fraction = ...
    estimated_spmv_total / baseline_time;


%% ========================================================================
% 15. Approximate OUTER FGMRES orthogonalisation model
%
% For q outer basis vectors, a classical full orthogonalisation requires
% approximately
%
%       q(q-1)/2
%
% vector-pair orthogonalisation interactions.
%
% Each interaction is approximated here by:
%
%       one dot product
%       one AXPY update
%
% plus one norm and one scaling operation per outer iteration.
%
% This is a MODEL, not a direct timing measurement.
% ========================================================================

orth_pair_count = ...
    outer_iterations * ...
    (outer_iterations - 1) / 2;


estimated_outer_orth_time = ...
    orth_pair_count * ...
    (dot_time + axpy_time) ...
    + ...
    outer_iterations * ...
    (norm_time + scale_time);


estimated_outer_orth_fraction = ...
    estimated_outer_orth_time / baseline_time;


fprintf('\n');
fprintf('---------------- MICROBENCHMARK RESULTS ----------------\n');

fprintf('Mean/typical SpMV cost       : %.8e s\n', ...
    spmv_time);

fprintf('Mean/typical dot cost        : %.8e s\n', ...
    dot_time);

fprintf('Mean/typical AXPY cost       : %.8e s\n', ...
    axpy_time);

fprintf('Mean/typical norm cost       : %.8e s\n', ...
    norm_time);

fprintf('Mean/typical scaling cost    : %.8e s\n', ...
    scale_time);


fprintf('\n');

fprintf('Estimated cumulative SpMV time : %.6f s\n', ...
    estimated_spmv_total);

fprintf('Estimated SpMV / baseline       : %.2f %%\n', ...
    100*estimated_spmv_fraction);


fprintf('\n');

fprintf('Outer orthogonalisation pairs   : %.0f\n', ...
    orth_pair_count);

fprintf('Modelled outer orth. time       : %.6f s\n', ...
    estimated_outer_orth_time);

fprintf('Modelled outer orth. / baseline : %.2f %%\n', ...
    100*estimated_outer_orth_fraction);

fprintf('---------------------------------------------------------\n\n');


%% ========================================================================
% 16. RUN 2:
%     MATLAB PROFILER
%
% The purpose of this run is NOT accurate end-to-end timing.
%
% It identifies:
%       - expensive functions
%       - expensive lines inside fastgmres.m
%
% Profiling introduces overhead, therefore this run MUST NOT replace
% baseline_time.
% ========================================================================

fprintf('============================================================\n');
fprintf(' RUN 2: MATLAB PROFILER\n');
fprintf('============================================================\n\n');


fprintf(['This run may take longer than RUN 1 because MATLAB is ' ...
         'collecting profiling information.\n\n']);


profile clear;

profile on;


TOTAL_SPMV_COUNT = 0;


% Reset randomized state so the solver receives the same sketch sequence
rng('default');


warning_state = warning;

warning('off','all');


tic;


[x_profile, ...
 rv_profile, ...
 rt_profile, ...
 spmv_profile] = ...
    fastgmres( ...
        A, ...
        b, ...
        tol);


profile_wall_time = toc;


warning(warning_state);


profile off;


prof_info = ...
    profile('info');


profile_relres = ...
    norm(b-A*x_profile)/norm_b;


fprintf('Profiler run completed.\n');

fprintf('Profiled wall-clock time : %.6f s\n', ...
    profile_wall_time);

fprintf('Profiled true residual   : %.8e\n\n', ...
    profile_relres);


fprintf(['IMPORTANT: the profiled wall-clock time above is NOT used ' ...
         'as the CPU baseline.\n']);

fprintf(['The accurate baseline remains RUN 1 = %.6f s.\n\n'], ...
    baseline_time);


%% ========================================================================
% 17. Print TOP profiler functions
% ========================================================================

fprintf('============================================================\n');
fprintf(' TOP FUNCTIONS REPORTED BY MATLAB PROFILER\n');
fprintf('============================================================\n\n');


FT = prof_info.FunctionTable;


if isempty(FT)

    fprintf('No profiler function information was returned.\n');

else

    function_times = ...
        [FT.TotalTime];


    [~,function_order] = ...
        sort( ...
            function_times, ...
            'descend');


    num_show_functions = ...
        min(20,length(function_order));


    fprintf('%-4s %-12s %-10s %s\n', ...
        'Rank', ...
        'Time (s)', ...
        'Calls', ...
        'Function');


    fprintf('%s\n', ...
        repmat('-',1,78));


    for kk = 1:num_show_functions

        idx = function_order(kk);


        fname = FT(idx).FunctionName;


        if isempty(fname)

            fname = '(unnamed)';

        end


        fprintf('%-4d %-12.6f %-10d %s\n', ...
            kk, ...
            FT(idx).TotalTime, ...
            FT(idx).NumCalls, ...
            fname);

    end

end


fprintf('\n');


%% ========================================================================
% 18. Print most expensive lines inside fastgmres.m
% ========================================================================

fprintf('============================================================\n');
fprintf(' HOT LINES INSIDE FASTGMRES.M\n');
fprintf('============================================================\n\n');


hot_line_data = [];


for k = 1:length(FT)

    file_name_lower = ...
        lower(FT(k).FileName);


    function_name_lower = ...
        lower(FT(k).FunctionName);


    is_fastgmres_entry = ...
        contains(file_name_lower,'fastgmres.m') ...
        || ...
        contains(function_name_lower,'fastgmres');


    if is_fastgmres_entry

        E = FT(k).ExecutedLines;


        if ~isempty(E) && size(E,2) >= 3

            % MATLAB ExecutedLines:
            %
            % column 1 = source line
            % column 2 = number of calls
            % column 3 = measured time

            nE = size(E,1);


            for jj = 1:nE

                hot_line_data = [
                    hot_line_data
                    E(jj,1), E(jj,2), E(jj,3), k
                ]; %#ok<AGROW>

            end

        end

    end

end


if isempty(hot_line_data)

    fprintf(['No line-level fastgmres information was available ' ...
             'from the profiler.\n\n']);

else

    [~,line_order] = ...
        sort( ...
            hot_line_data(:,3), ...
            'descend');


    num_show_lines = ...
        min(25,length(line_order));


    fprintf('%-4s %-8s %-10s %-12s %s\n', ...
        'Rank', ...
        'Line', ...
        'Calls', ...
        'Time (s)', ...
        'Function');


    fprintf('%s\n', ...
        repmat('-',1,90));


    for kk = 1:num_show_lines

        row = ...
            hot_line_data(line_order(kk),:);


        source_line = row(1);

        calls = row(2);

        line_time = row(3);

        function_idx = row(4);


        fprintf('%-4d %-8d %-10d %-12.6f %s\n', ...
            kk, ...
            source_line, ...
            calls, ...
            line_time, ...
            FT(function_idx).FunctionName);

    end

end


fprintf('\n');


%% ========================================================================
% 19. FINAL SUMMARY FOR CHAPTER 6 ANALYSIS
%
% THIS IS THE IMPORTANT BLOCK TO SCREENSHOT.
% ========================================================================

fprintf('\n');
fprintf('################################################################\n');
fprintf('#                                                              #\n');
fprintf('#        FINAL DATA FOR GPU SPEEDUP ESTIMATION                 #\n');
fprintf('#                                                              #\n');
fprintf('################################################################\n\n');


fprintf('A. ORIGINAL CPU FASTGMRES SOLVE\n');

fprintf('------------------------------------------------------------\n');

fprintf('Baseline CPU runtime             = %.6f s\n', ...
    baseline_time);

fprintf('True relative residual           = %.8e\n', ...
    true_relres);

fprintf('Outer iterations                 = %d\n', ...
    outer_iterations);

fprintf('Total SpMVs                      = %.0f\n', ...
    total_spmv);

fprintf('Estimated inner SpMVs            = %.0f\n', ...
    inner_spmv);

fprintf('Average inner steps per outer    = %.3f\n', ...
    average_inner_steps);


fprintf('\n');


fprintf('B. REPRESENTATIVE CPU KERNEL COSTS\n');

fprintf('------------------------------------------------------------\n');

fprintf('One SpMV                         = %.8e s\n', ...
    spmv_time);

fprintf('One dot product                  = %.8e s\n', ...
    dot_time);

fprintf('One AXPY                         = %.8e s\n', ...
    axpy_time);

fprintf('One vector norm                  = %.8e s\n', ...
    norm_time);

fprintf('One vector scaling               = %.8e s\n', ...
    scale_time);


fprintf('\n');


fprintf('C. MICROBENCHMARK-BASED WORKLOAD ESTIMATES\n');

fprintf('------------------------------------------------------------\n');

fprintf('Estimated cumulative SpMV time   = %.6f s\n', ...
    estimated_spmv_total);

fprintf('Estimated SpMV runtime fraction  = %.2f %%\n', ...
    100*estimated_spmv_fraction);

fprintf('Outer orthogonalisation pairs    = %.0f\n', ...
    orth_pair_count);

fprintf('Modelled outer orth. time        = %.6f s\n', ...
    estimated_outer_orth_time);

fprintf('Modelled outer orth. fraction    = %.2f %%\n', ...
    100*estimated_outer_orth_fraction);


fprintf('\n');


fprintf('D. PROFILER INFORMATION\n');

fprintf('------------------------------------------------------------\n');

fprintf('Profiler wall time               = %.6f s\n', ...
    profile_wall_time);

fprintf('(DO NOT use this as baseline runtime.)\n');


fprintf('\n');

fprintf(['Please also capture the TOP FUNCTIONS and HOT LINES ' ...
         'tables printed immediately above.\n']);


fprintf('\n');

fprintf('################################################################\n');
fprintf('# END OF PROFILING EXPERIMENT                                  #\n');
fprintf('################################################################\n');


%% ========================================================================
% LOCAL FUNCTION 1
% Benchmark sparse matrix-vector multiplication.
%
% Returns median block-average cost per operation.
% ========================================================================

function t = benchmark_spmv(A,x,nBlocks,nRepeat)

    block_times = zeros(nBlocks,1);

    y = A*x; %#ok<NASGU>


    for b = 1:nBlocks

        tic;

        for k = 1:nRepeat

            y = A*x; %#ok<NASGU>

        end

        block_times(b) = toc/nRepeat;

    end


    t = median(block_times);

end


%% ========================================================================
% LOCAL FUNCTION 2
% Benchmark dot product.
% ========================================================================

function t = benchmark_dot(x,y,nBlocks,nRepeat)

    block_times = zeros(nBlocks,1);

    s = dot(x,y); %#ok<NASGU>


    for b = 1:nBlocks

        tic;

        for k = 1:nRepeat

            s = dot(x,y); %#ok<NASGU>

        end

        block_times(b) = toc/nRepeat;

    end


    t = median(block_times);

end


%% ========================================================================
% LOCAL FUNCTION 3
% Benchmark AXPY-style vector operation:
%
%       z = x + alpha*y
% ========================================================================

function t = benchmark_axpy(x,y,alpha,nBlocks,nRepeat)

    block_times = zeros(nBlocks,1);

    z = x + alpha*y; %#ok<NASGU>


    for b = 1:nBlocks

        tic;

        for k = 1:nRepeat

            z = x + alpha*y; %#ok<NASGU>

        end

        block_times(b) = toc/nRepeat;

    end


    t = median(block_times);

end


%% ========================================================================
% LOCAL FUNCTION 4
% Benchmark vector 2-norm.
% ========================================================================

function t = benchmark_norm(x,nBlocks,nRepeat)

    block_times = zeros(nBlocks,1);

    s = norm(x); %#ok<NASGU>


    for b = 1:nBlocks

        tic;

        for k = 1:nRepeat

            s = norm(x); %#ok<NASGU>

        end

        block_times(b) = toc/nRepeat;

    end


    t = median(block_times);

end


%% ========================================================================
% LOCAL FUNCTION 5
% Benchmark vector scaling:
%
%       z = alpha*x
% ========================================================================

function t = benchmark_scale(x,alpha,nBlocks,nRepeat)

    block_times = zeros(nBlocks,1);

    z = alpha*x; %#ok<NASGU>


    for b = 1:nBlocks

        tic;

        for k = 1:nRepeat

            z = alpha*x; %#ok<NASGU>

        end

        block_times(b) = toc/nRepeat;

    end


    t = median(block_times);

end