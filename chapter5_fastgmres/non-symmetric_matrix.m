% =========================================================================
% Script: Load FluAtm Non-Symmetric CFD Matrix and Plot Sparsity Pattern
% Description: Reads binary files exported from the C/CUDA bridge,
%              reconstructs the sparse matrix A, and exports a 
%              high-resolution plot using custom supervisor configurations.
% =========================================================================
clear; clc; close all;

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

fprintf('🚀 [MATLAB] Loading FluAtm Non-Symmetric Matrix...\n');

%% 1. Read Row Pointer Array
fid = fopen('fluatm_matrix_row.dat', 'rb');
if fid == -1, error('Failed to open fluatm_matrix_row.dat'); end
n_row_ptr = fread(fid, 1, 'int32');
row_ptr = fread(fid, n_row_ptr, 'int32');
fclose(fid);
nrows = n_row_ptr - 1; 

%% 2. Read Column Indices Array
fid = fopen('fluatm_matrix_col.dat', 'rb');
if fid == -1, error('Failed to open fluatm_matrix_col.dat'); end
nnz = fread(fid, 1, 'int32');
col_idx = fread(fid, nnz, 'int32');
fclose(fid);

%% 3. Read Non-Zero Numerical Values Array
fid = fopen('fluatm_matrix_val.dat', 'rb');
if fid == -1, error('Failed to open fluatm_matrix_val.dat'); end
nnz_val = fread(fid, 1, 'int32');
values = fread(fid, nnz_val, 'double');
fclose(fid);

%% 4. Read Right-Hand Side (RHS) Vector
fid = fopen('fluatm_rhs.dat', 'rb');
if fid == -1, error('Failed to open fluatm_rhs.dat'); end
x_len = fread(fid, 1, 'int32');
b = fread(fid, x_len, 'double');
fclose(fid);

fprintf('✅ Data loaded successfully. Dimension: %d x %d, NNZ: %d\n', nrows, nrows, nnz);

%% 5. Convert CSR Format to COO Coordinates
fprintf('🔧 Converting CSR format to MATLAB sparse coordinate format...\n');
row_counts = diff(row_ptr);
row_idx = repelem((1:nrows)', row_counts);
A = sparse(row_idx, col_idx, values, nrows, nrows);

fprintf('======================================================\n');
fprintf('[Matrix Properties]\n');
fprintf('- Matrix Size:        %d x %d\n', nrows, nrows);
fprintf('- Non-Zero Elements:  %d\n', nnz);
fprintf('- Is Symmetric?       %s\n', mat2str(issymmetric(A)));
fprintf('======================================================\n');

%% 6. Sparsity Pattern Plotting with Supervisor Configurations
fig = figure('Name', 'FluAtm Matrix Sparsity Pattern', 'Color', 'w'); 

spy(A, '.', 0.2); 

% Set the color of the markers to an academic dark red for non-symmetric systems
h = get(gca, 'Children');
set(h, 'Color', [0.6, 0.0, 0.0]); 

% Adjust the axes plotting box to create comfortable margins at the top and left
ax = gca;
ax.Position(4) = ax.Position(4) - 0.05; 

% Force tick marks, tick labels, and boundary lines of the axes to be solid black
set(gca, 'XColor', 'k', 'YColor', 'k');

% Force the title to float safely above the scientific notation (x 10^4) using normalized coordinates
title('Sparsity Pattern of the FluAtm Momentum Matrix', ...
      'Units', 'normalized', ...
      'Position', [0.5, 1.12, 0], ...
      'FontWeight', 'bold', ...
      'Color', 'k');

xlabel(sprintf('Columns (N = %d)', nrows), 'Color', 'k');
ylabel(sprintf('Rows (N = %d)', nrows), 'Color', 'k');

% Force a quick render update and brief pause to ensure graphics queue is stable
drawnow;
pause(0.2);

% Export using custom template function
mypdf(fig, 'FluAtm_Matrix_Sparsity'); 
disp('[MATLAB] Matrix successfully processed and exported!');


function mypdf(fig, fname, r, s, pngonly)
% fig: figure handle to print
% r: height/width ratio, default  r = .66
% s: scaling of font size, default s = 1.5
% pngonly: only output png, default = 0

if nargin < 3, r = .66; end
if nargin < 4, s = 1.5; end
if nargin < 5, pngonly = 0; end

% Format paper size on the specified figure handle (for .fig compatibility)
if isgraphics(fig)
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperSize', s * [13, r * 13]);
    set(fig, 'PaperPosition', s * [0, 0, 13, r * 13]);
end

% Robust error-shielding blocks (try-catch) combined with active graphics handle checks.
% In headless Linux server environments, the first export call might silently close
% the figure's active graphics peer, making subsequent print commands fail.
% This design ensures that the script runs successfully without crashing.

if ~pngonly && isgraphics(fig)
    try
        exportgraphics(fig, [fname '.pdf'], 'ContentType', 'vector');
    catch
        warning('Unable to write vector PDF format in current headless environment.');
    end
end

if ~pngonly && isgraphics(fig)
    try
        exportgraphics(fig, [fname '.eps'], 'ContentType', 'vector');
    catch
        warning('Unable to write vector EPS format in current headless environment.');
    end
end

if isgraphics(fig)
    try
        exportgraphics(fig, [fname '.png'], 'Resolution', 200);
    catch
        warning('Unable to write raster PNG format in current headless environment.');
    end
end

if ~pngonly && isgraphics(fig)
    try
        saveas(fig, [fname '.fig']);
    catch
        warning('Unable to write native FIG format.');
    end
end
end