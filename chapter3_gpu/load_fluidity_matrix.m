% =========================================================================
% Script: Plot_Fluidity_Pressure_Matrix_PDF.m
%
% Purpose:
%   Read the Fluidity pressure matrix exported in CSR format and create
%   a publication-quality sparsity-pattern PDF.
%
% Key idea:
%   Instead of asking MATLAB to export millions of individual SPY markers,
%   the nonzero pattern is mapped to a high-resolution occupancy image.
%
% Output:
%   Fluidity_Pressure_Matrix_Sparsity.pdf
%   Fluidity_Pressure_Matrix_Sparsity_preview.png
% =========================================================================

clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf(' FLUIDITY PRESSURE MATRIX SPARSITY PLOT\n');
fprintf('============================================================\n\n');


%% ========================================================================
% 1. Input files
% ========================================================================

row_file = 'fluidity_matrix_row.dat';
col_file = 'fluidity_matrix_col.dat';
val_file = 'fluidity_matrix_val.dat';
rhs_file = 'fluidity_rhs.dat';

files = {
    row_file
    col_file
    val_file
    rhs_file
};

fprintf('Checking input files...\n');

for k = 1:length(files)

    if ~isfile(files{k})
        error('Required input file "%s" was not found.', files{k});
    end

    fprintf('  %-30s : OK\n', files{k});

end

fprintf('\n');


%% ========================================================================
% 2. Read CSR row pointer
% ========================================================================

fprintf('[1/4] Reading row pointers...\n');

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

fprintf('      Matrix dimension N = %d\n',nrows);


%% ========================================================================
% 3. Read CSR column indices
% ========================================================================

fprintf('[2/4] Reading column indices...\n');

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

fprintf('      Column entries = %d\n',nnz_col);


%% ========================================================================
% 4. Read numerical values
% ========================================================================

fprintf('[3/4] Reading numerical values...\n');

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

fprintf('      Numerical values = %d\n',nnz_val);


%% ========================================================================
% 5. Read RHS
% ========================================================================

fprintf('[4/4] Reading RHS vector...\n');

fid = fopen(rhs_file,'rb');

if fid == -1
    error('Failed to open %s.',rhs_file);
end

rhs_len = fread(fid,1,'int32');
b = fread(fid,rhs_len,'double');

fclose(fid);

fprintf('      RHS length = %d\n\n',rhs_len);


%% ========================================================================
% 6. Consistency checks
% ========================================================================

if nnz_col ~= nnz_val

    error( ...
        'Column-index count (%d) differs from value count (%d).', ...
        nnz_col,nnz_val);

end

if rhs_len ~= nrows

    warning( ...
        'RHS length (%d) differs from matrix dimension (%d).', ...
        rhs_len,nrows);

end


%% ========================================================================
% 7. Detect column-index convention
% ========================================================================

fprintf('Checking column-index convention...\n');

if min(col_idx) == 0

    fprintf('  Detected 0-based indexing.\n');
    fprintf('  Converting to MATLAB 1-based indexing.\n');

    col_idx = col_idx + 1;

elseif min(col_idx) >= 1

    fprintf('  Column indices already appear to be 1-based.\n');

else

    error('Negative column index detected.');

end

if max(col_idx) > nrows
    error('Column index exceeds matrix dimension.');
end

fprintf('\n');


%% ========================================================================
% 8. Convert CSR row pointers into explicit row coordinates
% ========================================================================

fprintf('Converting CSR structure to row/column coordinates...\n');

row_counts = diff(row_ptr);

if any(row_counts < 0)
    error('CSR row pointers are not monotonically increasing.');
end

if sum(double(row_counts)) ~= nnz_val

    error( ...
        'CSR structure predicts %d entries but value file contains %d.', ...
        sum(double(row_counts)),nnz_val);

end

row_idx = repelem( ...
    (1:nrows)', ...
    double(row_counts));

fprintf('Coordinate conversion completed.\n\n');


%% ========================================================================
% 9. Matrix information
% ========================================================================

fprintf('============================================================\n');
fprintf(' MATRIX INFORMATION\n');
fprintf('============================================================\n');

fprintf('Dimension : %d x %d\n', ...
    nrows,nrows);

fprintf('NNZ       : %d\n', ...
    nnz_val);

fprintf('Density   : %.6e\n', ...
    double(nnz_val)/(double(nrows)*double(nrows)));

fprintf('============================================================\n\n');


%% ========================================================================
% 10. Construct high-resolution occupancy image
%
% IMPORTANT:
%
% The original matrix is 95187 x 95187, much larger than any printed
% dissertation figure.
%
% We map the sparsity pattern onto a 2500 x 2500 image. A pixel is dark
% blue if one or more matrix nonzeros fall inside the corresponding region.
%
% This preserves the visible sparsity structure while avoiding millions of
% MATLAB graphical marker objects.
% ========================================================================

image_resolution = 2500;

fprintf('Generating %d x %d sparsity image...\n', ...
    image_resolution,image_resolution);


% -------------------------------------------------------------------------
% Map original matrix coordinates to image pixels
% -------------------------------------------------------------------------

pixel_row = ...
    floor( ...
        (double(row_idx)-1) ...
        * image_resolution ...
        / double(nrows)) ...
    + 1;


pixel_col = ...
    floor( ...
        (double(col_idx)-1) ...
        * image_resolution ...
        / double(nrows)) ...
    + 1;


% Numerical protection
pixel_row = min(max(pixel_row,1),image_resolution);
pixel_col = min(max(pixel_col,1),image_resolution);


% -------------------------------------------------------------------------
% Build occupancy map
%
% sparse() is used here first because it is efficient for millions of
% coordinate entries.
% -------------------------------------------------------------------------

pattern_sparse = sparse( ...
    pixel_row, ...
    pixel_col, ...
    1, ...
    image_resolution, ...
    image_resolution);


% Convert to logical occupancy
pattern = full(pattern_sparse > 0);

clear pattern_sparse
clear pixel_row pixel_col
clear row_idx col_idx values b

fprintf('Sparsity image generated successfully.\n\n');


%% ========================================================================
% 11. Create publication-quality figure
% ========================================================================

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 7.2 6.6]);


ax = axes(fig);


%% ========================================================================
% 12. Display sparsity image
% ========================================================================

imagesc( ...
    ax, ...
    [1 nrows], ...
    [1 nrows], ...
    pattern);


axis(ax,'square');


% Matrix convention:
% row 1 should appear at the top
set(ax,'YDir','reverse');


%% ========================================================================
% 13. Colour map
%
% 0 = white
% 1 = dark academic blue
% ========================================================================

colormap( ...
    ax, ...
    [
        1.00 1.00 1.00
        0.00 0.20 0.40
    ]);


caxis(ax,[0 1]);


%% ========================================================================
% 14. Axes formatting
% ========================================================================

ax.Color = 'w';

ax.XColor = 'k';
ax.YColor = 'k';

ax.FontName = 'Times New Roman';
ax.FontSize = 11;

ax.LineWidth = 0.9;

grid(ax,'off');
box(ax,'on');


xlim(ax,[1 nrows]);
ylim(ax,[1 nrows]);


xlabel( ...
    ax, ...
    sprintf('Column index (N = %d)',nrows), ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'Color','k');


ylabel( ...
    ax, ...
    sprintf('Row index (N = %d)',nrows), ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'Color','k');


title( ...
    ax, ...
    'Fluidity Pressure Matrix Sparsity Pattern', ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'FontWeight','normal', ...
    'Color','k');


%% ========================================================================
% 15. Adjust spacing
% ========================================================================

ax.Position = ...
    [0.14 0.14 0.78 0.76];


drawnow;


%% ========================================================================
% 16. Export preview PNG
% ========================================================================

preview_file = ...
    'Fluidity_Pressure_Matrix_Sparsity_preview.png';


fprintf('Exporting preview PNG...\n');


exportgraphics( ...
    fig, ...
    preview_file, ...
    'Resolution', ...
    300, ...
    'BackgroundColor', ...
    'white');


fprintf('Preview PNG completed.\n');


%% ========================================================================
% 17. Export PDF
%
% Because the plotted object is now one image rather than millions of
% individual SPY points, this export is much lighter and more reliable.
% ========================================================================

pdf_file = ...
    'Fluidity_Pressure_Matrix_Sparsity.pdf';


fprintf('Exporting PDF...\n');


try

    exportgraphics( ...
        fig, ...
        pdf_file, ...
        'ContentType', ...
        'image', ...
        'Resolution', ...
        400, ...
        'BackgroundColor', ...
        'white');

    fprintf('PDF export completed using exportgraphics.\n');


catch ME

    fprintf('\n');
    fprintf('exportgraphics PDF export failed:\n');
    fprintf('%s\n\n',ME.message);

    fprintf('Trying MATLAB print() fallback...\n');


    % ---------------------------------------------------------------------
    % Fallback export path
    % ---------------------------------------------------------------------

    set(fig,'PaperPositionMode','auto');

    set(fig,'Renderer','opengl');


    print( ...
        fig, ...
        pdf_file, ...
        '-dpdf', ...
        '-opengl', ...
        '-r400');


    fprintf('PDF export completed using print().\n');

end


%% ========================================================================
% 18. Final message
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' EXPORT COMPLETED\n');
fprintf('============================================================\n\n');

fprintf('Dissertation PDF:\n');
fprintf('  %s\n\n',pdf_file);

fprintf('Local preview:\n');
fprintf('  %s\n\n',preview_file);

fprintf('Upload the PDF directly to Overleaf.\n');

fprintf('\n============================================================\n');
fprintf(' DONE\n');
fprintf('============================================================\n');