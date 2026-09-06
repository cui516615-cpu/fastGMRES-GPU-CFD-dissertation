% =========================================================================
% Script: Run_fastGMRES_Comparison.m
% Description: Matched time-budget benchmarking with true SpMV tracking.
% =========================================================================

clear;
clc;
close all;

% Set default styling configurations for publication-quality figures
randn('state', 0);
rand('state', 0);

height = 2/3;

set(0, ...
    'defaultfigureposition', [180 100 800 800*height], ...
    'defaultaxeslinewidth', 1, ...
    'defaultaxesfontsize', 18, ...
    'defaultlinelinewidth', 2, ...
    'defaultpatchlinewidth', 1, ...
    'defaultlinemarkersize', 10, ...
    'defaulttextinterpreter', 'tex');


fprintf('🚀 [Phase 1] Ingesting binary matrix data...\n');


%% =========================================================================
% Read row pointers
% =========================================================================

fid = fopen('fluatm_matrix_row.dat', 'rb');

if fid == -1
    error('Failed to open fluatm_matrix_row.dat');
end

n_row_ptr = fread(fid, 1, 'int32');
row_ptr = fread(fid, n_row_ptr, 'int32');

fclose(fid);

nrows = n_row_ptr - 1;


%% =========================================================================
% Read column indices
% =========================================================================

fid = fopen('fluatm_matrix_col.dat', 'rb');

if fid == -1
    error('Failed to open fluatm_matrix_col.dat');
end

nnz_col = fread(fid, 1, 'int32');
col_idx = fread(fid, nnz_col, 'int32');

fclose(fid);


%% =========================================================================
% Read matrix values
% =========================================================================

fid = fopen('fluatm_matrix_val.dat', 'rb');

if fid == -1
    error('Failed to open fluatm_matrix_val.dat');
end

nnz_val = fread(fid, 1, 'int32');
values = fread(fid, nnz_val, 'double');

fclose(fid);


%% =========================================================================
% Read RHS
% =========================================================================

fid = fopen('fluatm_rhs.dat', 'rb');

if fid == -1
    error('Failed to open fluatm_rhs.dat');
end

x_len = fread(fid, 1, 'int32');
b = fread(fid, x_len, 'double');

fclose(fid);


%% =========================================================================
% Assemble sparse matrix
% =========================================================================

A = sparse( ...
    repelem((1:nrows)', diff(row_ptr)), ...
    col_idx, ...
    values, ...
    nrows, ...
    nrows);


norm_b = norm(b);


fprintf('✅ System loaded successfully.\n');

fprintf('======================================================\n');
fprintf('🚀 [Phase 2] Executing Solver Benchmarks...\n');


%% =========================================================================
% Solver settings
% =========================================================================

tol = 1e-7;

restart = 50;

max_outer = 120;

max_iter_bicg = 6000;


%% =========================================================================
% 1. Standard GMRES(50)
% =========================================================================

fprintf('Running Standard GMRES...\n');

tic;

[x_gmres, ...
 flag_gmres, ...
 relres_gmres, ...
 iter_gmres, ...
 rv_gmres] = ...
    gmres( ...
        A, ...
        b, ...
        restart, ...
        tol, ...
        max_outer);

time_gmres = toc;


%% =========================================================================
% 2. BiCGSTAB
% =========================================================================

fprintf('Running BiCGSTAB...\n');

tic;

[x_bicg, ...
 flag_bicg, ...
 relres_bicg, ...
 iter_bicg, ...
 rv_bicg] = ...
    bicgstab( ...
        A, ...
        b, ...
        tol, ...
        max_iter_bicg);

time_bicg = toc;


%% =========================================================================
% 3. fastGMRES
% =========================================================================

global TOTAL_SPMV_COUNT;

TOTAL_SPMV_COUNT = 0;


fprintf('Running fastGMRES with true SpMV tracking...\n');


warning('off', 'all');

tic;

[x_fast, ...
 rv_fast, ...
 rt_fast, ...
 spmv_history] = ...
    fastgmres( ...
        A, ...
        b, ...
        tol);

time_fast = toc;

warning('on', 'all');


relres_fast = ...
    norm(b - A*x_fast) / norm_b;


%% =========================================================================
% Phase 3: Convergence Plot
%
% IMPORTANT:
% The plotting logic below is intentionally kept the same as the original.
% Only the background colours and PDF export settings are changed.
% =========================================================================

fprintf('======================================================\n');

fprintf('🎨 Exporting convergence curves via mypdf...\n');


%% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------

fig = figure( ...
    'Name', ...
    'Convergence History Comparison', ...
    'Color', ...
    'w');


%% -------------------------------------------------------------------------
% GMRES(50)
%
% Original plotting rule:
% each inner iteration corresponds to one SpMV.
% -------------------------------------------------------------------------

semilogy( ...
    0:length(rv_gmres)-1, ...
    rv_gmres / norm_b, ...
    'b-', ...
    'LineWidth', ...
    2.5);

hold on;


%% -------------------------------------------------------------------------
% BiCGSTAB
%
% Original plotting rule:
% rv_bicg residual history is used directly against active SpMV count.
% -------------------------------------------------------------------------

semilogy( ...
    0:length(rv_bicg)-1, ...
    rv_bicg / norm_b, ...
    'r--', ...
    'LineWidth', ...
    2);


%% -------------------------------------------------------------------------
% fastGMRES
%
% Original plotting rule:
% outer residuals are plotted at exact cumulative SpMV coordinates.
% -------------------------------------------------------------------------

semilogy( ...
    spmv_history, ...
    rv_fast / norm_b, ...
    'Color', ...
    [0.0, 0.6, 0.3], ...
    'LineStyle', ...
    '-.', ...
    'LineWidth', ...
    3.0);


%% =========================================================================
% ORIGINAL axes formatting
% =========================================================================

set( ...
    gca, ...
    'XColor', ...
    'k', ...
    'YColor', ...
    'k');


ax = gca;


% -------------------------------------------------------------------------
% ONLY NEW CHANGE:
% force the plotting region to WHITE.
% -------------------------------------------------------------------------

ax.Color = 'w';


% Preserve logarithmic y-axis used by semilogy
ax.YScale = 'log';


% -------------------------------------------------------------------------
% Original spacing
% -------------------------------------------------------------------------

ax.Position(4) = ...
    ax.Position(4) - 0.05;


%% =========================================================================
% Original grid settings
% =========================================================================

grid on;

set(gca, ...
    'YMinorGrid', ...
    'on');


%% =========================================================================
% Original labels
% =========================================================================

xlabel( ...
    'Number of Matrix-Vector Products (SpMVs)', ...
    'Color', ...
    'k');


ylabel( ...
    'Relative Residual Norm: ||b - Ax|| / ||b||', ...
    'Color', ...
    'k');


%% =========================================================================
% Original title
% =========================================================================

title( ...
    'Convergence Comparison on FluAtm Momentum Matrix', ...
    'Units', ...
    'normalized', ...
    'Position', ...
    [0.5, 1.12, 0], ...
    'FontWeight', ...
    'bold', ...
    'Color', ...
    'k');


%% =========================================================================
% Original legend
%
% ONLY background colour is changed.
% =========================================================================

lgd = legend( ...
    'Standard GMRES (restart=50)', ...
    'BiCGSTAB', ...
    'fastGMRES (Proposed)', ...
    'Location', ...
    'northeast');


% -------------------------------------------------------------------------
% ONLY NEW CHANGES:
% white legend background and black text.
% -------------------------------------------------------------------------

lgd.Color = 'w';

lgd.TextColor = 'k';

lgd.EdgeColor = 'k';


%% =========================================================================
% Force whole figure background to white
% =========================================================================

fig.Color = 'w';


drawnow;

pause(0.2);


%% =========================================================================
% Export
% =========================================================================

mypdf( ...
    fig, ...
    'FluAtm_Convergence_Comparison');


fprintf( ...
    '✅ Benchmarking complete. Comparison plots exported successfully.\n');


%% =========================================================================
% mypdf
%
% Same original function.
% Only BackgroundColor='white' has been added to exportgraphics.
% =========================================================================

function mypdf(fig, fname, r, s, pngonly)

    if nargin < 3
        r = .66;
    end

    if nargin < 4
        s = 1.5;
    end

    if nargin < 5
        pngonly = 0;
    end


    if isgraphics(fig)

        set( ...
            fig, ...
            'PaperPositionMode', ...
            'auto');


        set( ...
            fig, ...
            'PaperSize', ...
            s * [13, r * 13]);


        set( ...
            fig, ...
            'PaperPosition', ...
            s * [0, 0, 13, r * 13]);


        % Force white figure background
        set( ...
            fig, ...
            'Color', ...
            'w');

    end


    %% --------------------------------------------------------------------
    % PDF
    % ---------------------------------------------------------------------

    if ~pngonly && isgraphics(fig)

        try

            exportgraphics( ...
                fig, ...
                [fname '.pdf'], ...
                'ContentType', ...
                'vector', ...
                'BackgroundColor', ...
                'white');

        catch ME

            fprintf( ...
                'PDF export error: %s\n', ...
                ME.message);

        end

    end


    %% --------------------------------------------------------------------
    % EPS
    % ---------------------------------------------------------------------

    if ~pngonly && isgraphics(fig)

        try

            exportgraphics( ...
                fig, ...
                [fname '.eps'], ...
                'ContentType', ...
                'vector', ...
                'BackgroundColor', ...
                'white');

        catch ME

            fprintf( ...
                'EPS export error: %s\n', ...
                ME.message);

        end

    end


    %% --------------------------------------------------------------------
    % PNG preview
    % ---------------------------------------------------------------------

    if isgraphics(fig)

        try

            exportgraphics( ...
                fig, ...
                [fname '.png'], ...
                'Resolution', ...
                200, ...
                'BackgroundColor', ...
                'white');

        catch ME

            fprintf( ...
                'PNG export error: %s\n', ...
                ME.message);

        end

    end


    %% --------------------------------------------------------------------
    % MATLAB FIG
    % ---------------------------------------------------------------------

    if ~pngonly && isgraphics(fig)

        try

            saveas( ...
                fig, ...
                [fname '.fig']);

        catch ME

            fprintf( ...
                'FIG export error: %s\n', ...
                ME.message);

        end

    end

end