clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf(' RE-PLOT CORRECTED PERFORMANCE PROFILES\n');
fprintf('============================================================\n\n');


%% ========================================================================
% 1. Load existing benchmark results
%
% IMPORTANT:
% This script does NOT rerun any solver.
% It only uses the already saved benchmark data.
% ========================================================================

result_file = ...
    'PublicMatrix_20_CASE_ILU0_results.mat';


if ~isfile(result_file)

    error( ...
        'Result file "%s" was not found.', ...
        result_file);

end


load(result_file, ...
    'performance_time', ...
    'Results');


fprintf('Loaded saved benchmark results successfully.\n\n');


%% ========================================================================
% 2. Solver names
% ========================================================================

methods = {
    'GMRES(50)'
    'BiCGSTAB'
    'fastGMRES'
};


%% ========================================================================
% 3. Define benchmark groups
%
% Rows:
%
% 1  NACA12
% 2  NACA12-p
% 3  ML_Laplace
% 4  ML_Laplace-p
% ...
%
% Therefore:
%
% odd rows  = without preconditioning
% even rows = ILU(0)
% ========================================================================

numCases = ...
    size(performance_time,1);


rows_no_prec = ...
    1:2:numCases;


rows_ilu = ...
    2:2:numCases;


fprintf('Total benchmark cases          : %d\n', ...
    numCases);

fprintf('Unpreconditioned cases         : %d\n', ...
    length(rows_no_prec));

fprintf('ILU(0)-preconditioned cases    : %d\n\n', ...
    length(rows_ilu));


%% ========================================================================
% 4. Print TRUE failure counts
% ========================================================================

fprintf('============================================================\n');
fprintf(' TRUE FAILURE COUNTS\n');
fprintf('============================================================\n\n');


fail_all = ...
    sum(~isfinite(performance_time),1);


fail_no_prec = ...
    sum(~isfinite( ...
        performance_time(rows_no_prec,:)),1);


fail_ilu = ...
    sum(~isfinite( ...
        performance_time(rows_ilu,:)),1);


fprintf('ALL 20 CASES\n');

fprintf('GMRES(50) : %d failures\n', ...
    fail_all(1));

fprintf('BiCGSTAB   : %d failures\n', ...
    fail_all(2));

fprintf('fastGMRES  : %d failures\n\n', ...
    fail_all(3));


fprintf('WITHOUT PRECONDITIONING\n');

fprintf('GMRES(50) : %d failures\n', ...
    fail_no_prec(1));

fprintf('BiCGSTAB   : %d failures\n', ...
    fail_no_prec(2));

fprintf('fastGMRES  : %d failures\n\n', ...
    fail_no_prec(3));


fprintf('WITH ILU(0)\n');

fprintf('GMRES(50) : %d failures\n', ...
    fail_ilu(1));

fprintf('BiCGSTAB   : %d failures\n', ...
    fail_ilu(2));

fprintf('fastGMRES  : %d failures\n\n', ...
    fail_ilu(3));


%% ========================================================================
% 5. PERFORMANCE PROFILE:
% ALL 20 CASES
% ========================================================================

fig1 = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [200 150 900 600]);


pprofile_corrected( ...
    performance_time, ...
    methods);


title( ...
    'Runtime Performance Profile: 20 Benchmark Cases', ...
    'FontWeight', ...
    'normal', ...
    'Color', ...
    'k');


exportgraphics( ...
    fig1, ...
    'Performance_Profile_20Cases_CORRECTED.pdf', ...
    'ContentType', ...
    'vector', ...
    'BackgroundColor', ...
    'white');


%% ========================================================================
% 6. PERFORMANCE PROFILE:
% WITHOUT PRECONDITIONING
% ========================================================================

fig2 = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [200 150 900 600]);


pprofile_corrected( ...
    performance_time(rows_no_prec,:), ...
    methods);


title( ...
    'Runtime Performance Profile: 10 Unpreconditioned Cases', ...
    'FontWeight', ...
    'normal', ...
    'Color', ...
    'k');


exportgraphics( ...
    fig2, ...
    'Performance_Profile_10_Unpreconditioned_CORRECTED.pdf', ...
    'ContentType', ...
    'vector', ...
    'BackgroundColor', ...
    'white');


%% ========================================================================
% 7. PERFORMANCE PROFILE:
% ILU(0)
% ========================================================================

fig3 = figure( ...
    'Color', ...
    'w', ...
    'Position', ...
    [200 150 900 600]);


pprofile_corrected( ...
    performance_time(rows_ilu,:), ...
    methods);


title( ...
    'Runtime Performance Profile: 10 ILU(0)-Preconditioned Cases', ...
    'FontWeight', ...
    'normal', ...
    'Color', ...
    'k');


exportgraphics( ...
    fig3, ...
    'Performance_Profile_10_ILU0_CORRECTED.pdf', ...
    'ContentType', ...
    'vector', ...
    'BackgroundColor', ...
    'white');


%% ========================================================================
% 8. Final message
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' CORRECTED FIGURES SAVED\n');
fprintf('============================================================\n\n');


fprintf( ...
    'Performance_Profile_20Cases_CORRECTED.pdf\n');

fprintf( ...
    'Performance_Profile_10_Unpreconditioned_CORRECTED.pdf\n');

fprintf( ...
    'Performance_Profile_10_ILU0_CORRECTED.pdf\n');


fprintf('\nDone.\n');



%% ========================================================================
% LOCAL FUNCTION
%
% Correct runtime performance profile
%
% IMPORTANT DIFFERENCE FROM THE OLD VERSION:
%
% ALL benchmark problems remain in the denominator.
%
% If all solvers fail on one problem:
%
%       [Inf Inf Inf]
%
% this problem is NOT removed.
%
% Instead, it contributes zero success to every solver while still
% remaining part of the total benchmark set.
%
% Therefore the final profile height also reflects solver robustness.
% ========================================================================

function pprofile_corrected(M,methods)


    %% --------------------------------------------------------------------
    % Number of ALL benchmark problems
    %
    % We DO NOT remove all-failure rows.
    % ---------------------------------------------------------------------

    nProblems = ...
        size(M,1);


    nMethods = ...
        size(M,2);


    fprintf('\nPerformance profile:\n');

    fprintf( ...
        '  Total benchmark cases = %d\n', ...
        nProblems);


    %% --------------------------------------------------------------------
    % Best successful runtime for each benchmark problem
    %
    % Example:
    %
    % [10 20 15]  -> best = 10
    %
    % [Inf 20 15] -> best = 15
    %
    % [Inf Inf Inf] -> best = Inf
    %
    % The final case remains in the denominator.
    % ---------------------------------------------------------------------

    best = ...
        min(M,[],2);


    %% --------------------------------------------------------------------
    % Compute finite performance ratios
    %
    % These are used only to choose a useful x-axis range.
    % ---------------------------------------------------------------------

    ratio = ...
        Inf(size(M));


    for j = 1:nMethods

        valid = ...
            isfinite(M(:,j)) & ...
            isfinite(best);

        ratio(valid,j) = ...
            M(valid,j)./best(valid);

    end


    finite_ratios = ...
        ratio(isfinite(ratio));


    %% --------------------------------------------------------------------
    % Choose x-axis range
    %
    % We want the right side of the plot to be large enough that all
    % successful solver cases appear in the profile.
    %
    % A minimum alpha_max = 4 is retained.
    % ---------------------------------------------------------------------

    if isempty(finite_ratios)

        alpha_max = 4;

    else

        alpha_max = ...
            max( ...
                4, ...
                1.05*max(finite_ratios));

    end


    % Round upward slightly for cleaner plotting
    alpha_max = ...
        ceil(alpha_max*10)/10;


    fprintf( ...
        '  Maximum performance ratio shown = %.1f\n', ...
        alpha_max);


    %% --------------------------------------------------------------------
    % Alpha values
    % ---------------------------------------------------------------------

    alpha = ...
        linspace(1,alpha_max,300);


    P = ...
        zeros(length(alpha),nMethods);


    %% --------------------------------------------------------------------
    % Construct performance profile
    %
    % KEY POINT:
    %
    % denominator = nProblems
    %
    % ALWAYS.
    %
    % Failed methods and all-failed problems therefore reduce the final
    % profile height.
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
            1.6, ...
            'MarkerIndices', ...
            1:30:length(alpha));

    end


    %% --------------------------------------------------------------------
    % White publication-style background
    % ---------------------------------------------------------------------

    ax = gca;

    ax.Color = ...
        'w';

    ax.XColor = ...
        'k';

    ax.YColor = ...
        'k';

    ax.GridColor = ...
        [0.7 0.7 0.7];

    ax.GridAlpha = ...
        0.35;

    ax.FontSize = ...
        11;


    grid on;

    box on;


    xlim( ...
        [1 alpha_max]);

    ylim( ...
        [0 1.05]);


    xlabel( ...
        'Performance ratio \alpha', ...
        'Color', ...
        'k');


    ylabel( ...
        'Fraction of test problems', ...
        'Color', ...
        'k');


    %% --------------------------------------------------------------------
    % TRUE failure counts
    %
    % Failures are counted over ALL benchmark cases.
    % ---------------------------------------------------------------------

    failures = ...
        sum(~isfinite(M),1);


    leg = ...
        cell(1,nMethods);


    for j = 1:nMethods

        if failures(j) == 0

            leg{j} = ...
                methods{j};

        elseif failures(j) == 1

            leg{j} = ...
                sprintf( ...
                    '%s (1 failure)', ...
                    methods{j});

        else

            leg{j} = ...
                sprintf( ...
                    '%s (%d failures)', ...
                    methods{j}, ...
                    failures(j));

        end

    end


    lgd = legend( ...
        leg, ...
        'Location', ...
        'southeast');


    lgd.Color = ...
        'w';

    lgd.TextColor = ...
        'k';

    lgd.EdgeColor = ...
        'k';


    %% --------------------------------------------------------------------
    % Print final profile heights
    %
    % Because alpha_max includes every finite performance ratio,
    % these should equal the solver success rates.
    % ---------------------------------------------------------------------

    fprintf('\n');

    fprintf( ...
        '  Final profile height GMRES(50) = %.3f\n', ...
        P(end,1));

    fprintf( ...
        '  Final profile height BiCGSTAB  = %.3f\n', ...
        P(end,2));

    fprintf( ...
        '  Final profile height fastGMRES = %.3f\n', ...
        P(end,3));


end