function [data_norm, nanpxs, data_train, data_test] = ProcessAndSplitDataAuditoryGngGreenDA(fnC_fn, save_fn, parameter_class)
% ProcessAndSplitDataAuditoryGngGreenDA
%
% Loads chunked GRAB-DA dF/F stacks for the auditory go/no-go task,
% removes NaN pixels, clips extreme values, applies a bypassed Lucy-
% Richardson preprocessing step (non-negativity enforcement only),
% reconstructs the stacks, optionally denoises them with PCA, normalizes
% the data, and finally splits the result into train and test tensors for
% downstream fpCNMF analysis.
%
% This version was modified to accommodate the data acquisition and
% train/test split scheme used for the auditory go/no-go task.
%
% INPUTS
%   fnC_fn          : path to a .mat file containing variable 'fnC'
%                     fnC must be an N x 2 cell array of file paths
%                     Column 1 = training chunks
%                     Column 2 = test chunks
%
%   save_fn         : output filename for saving processed results
%                     If empty, processed variables are returned but not saved
%
%   parameter_class : class name or constructor handle used to instantiate
%                     the parameter object, which is then passed through
%                     loadobj(...)
%
% OUTPUTS
%   data_norm  : normalized data matrix, size = nonNaN pixels x total frames
%   nanpxs     : indices of removed NaN pixels
%   data_train : training tensor,
%                size = nonNaN pixels x frames per chunk x # train chunks
%   data_test  : test tensor,
%                size = nonNaN pixels x frames per chunk x # test chunks
%
% NOTES
%   - Each input dF/F stack is expected to be X x Y x T
%   - conditionDffMat converts between stack format and frame-by-pixel format
%   - clipDataC_percentile and bypassLucric are local helper functions
%     defined at the end of this file and do not need to exist as separate
%     functions on the MATLAB path
%
% Original framework: Camden MacDowell - timeless
% Modified: Junchol Park, Feb 2025

load(fnC_fn, 'fnC'); 
% fnC = cellfun(@(a) ConvertBucketToWinPath(a), fnC, 'un', 0); % to convert the Spock path back to Mac Path

if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

gp = loadobj(feval(parameter_class)); 

num_chunks = size(fnC, 1); % number of train/test chunk pairs
assert(size(fnC, 2)==2);   % must contain exactly two columns: train and test
assert(sum(cellfun(@isempty, fnC(:)))==0); 

%% Load dF/F stacks and convert each stack to frame-by-pixel format
for rr = 1:size(fnC, 1)
    for cc = 1:size(fnC, 2)
        temp = load(fnC{rr, cc});
        if ~exist('opts', 'var')
            opts = temp.opts; 
        end      
        % record original stack dimensions
        xC{rr, cc} = size(temp.dff, 1); 
        yC{rr, cc} = size(temp.dff, 2);  
        zC{rr, cc} = size(temp.dff, 3);  
        
        % keep original stack for reference
        dataCrs_org{rr, cc} = temp.dff;
        
        % convert stack to [frames x non-NaN pixels]
        dataC{rr, cc} = [];
        nanpxsC{rr, cc} = [];
        [dataC{rr, cc}, nanpxsC{rr, cc}] = conditionDffMat(temp.dff); % frame-by-pixel matrix
        clear temp; 
        fprintf(sprintf("Completed loading dff or row#%d and col#%d\n", rr, cc)); 
    end
end

% Clip extreme values globally across all chunks to reduce the influence of
% rare artifactual outliers before subsequent processing
dataC = clipDataC_percentile(dataC, 0.2, 99.8); % Mask extreme values

%% Apply bypassed Lucy-Richardson preprocessing pixelwise
% This step preserves the non-negativity enforcement from lucric.m
% but skips the actual deconvolution. Each column is a pixel trace over time.
for rr = 1:size(fnC, 1)
    for cc = 1:size(fnC, 2)
        for px = 1:size(dataC{rr, cc},2)
            dataC{rr, cc}(:,px) = bypassLucric(dataC{rr, cc}(:,px), gp.d_gamma, gp.d_smooth, gp.d_kernel);
        end
        clear temp; 
        fprintf(sprintf("Completed loading dff or row#%d and col#%d\n", rr, cc)); 
    end
end

%% Dimension sanity checks
minFrN = min(cell2mat(zC(:)));  % minimum number of frames across all chunks
rowN = unique(cell2mat(xC(:))); % image height (expected scalar, e.g. 64)
assert(isscalar(rowN));         % all stacks must have same number of rows
colN = unique(cell2mat(yC(:))); % image width (expected scalar, e.g. 64)
assert(isscalar(colN));         % all stacks must have same number of columns
assert(isscalar(unique(cellfun(@nansum, nanpxsC(:))))); % NaN mask must match across chunks
nanpxs = nanpxsC{1}; 

%% Curtail to common frame count and reconstruct stacks
% Truncate all chunks to the same number of frames so that they can be
% stacked consistently across train/test splits
dataC = cellfun(@(a) a(1:minFrN, :), dataC, 'UniformOutput', false); 

% Convert each chunk back to image-stack format: X x Y x frames
dataC = cellfun(@(a, b) conditionDffMat(a, b, [], [rowN, colN, minFrN]), dataC, nanpxsC, 'UniformOutput', false); 

% Concatenate all reconstructed stacks along the frame dimension
data = cat(3, dataC{:}); % final size ~ [64 x 64 x total frames]
% isequaln(dataC{2, 1}, data(:,:,861:2*860)) % sanity check (true because MATLAB is column-major)

%% Optional PCA denoising
% Applied after reconstruction and before final frame-by-pixel conversion
if gp.w_pca_denoise
    data = DenoisePCA(data); % data: 64 x 64 x N total frames
end

% Convert back to [total frames x non-NaN pixels]
[data,~] = conditionDffMat(data); % data: total frames x non-NaN pixels

%% Normalize data
% Depending on gp.w_normalization_method, normalization is performed either
% pixelwise or globally.
fprintf('\n\tPerforming %s normalization to %d value', gp.w_normalization_method, gp.w_norm_val);
switch gp.w_normalization_method
    case 'pixelwise' % normalize each pixel independently by its own percentile
        data_norm = NaN(size(data));
        for px = 1:size(data,2)
            data_norm(:,px) = (data(:,px))/(prctile(data(:,px),gp.w_norm_val));
        end             
    case 'full' % normalize all data by a global percentile computed from positive values
        data_norm = data/prctile(data(data>eps), gp.w_norm_val);          
    case 'bounded'
        data_norm = (data)/(gp.w_norm_val(2)); % normalize by provided upper bound
    case 'none'
        data_norm = data;
    otherwise
        error('Unknown normalization method. Check general params')
end

%% Transpose to match fpCNMF convention
% fpCNMF operates rowwise, so output should be:
% [non-NaN pixels x total frames]
data_norm = data_norm'; 

%% Create train and test tensors
assert(isequal(size(data_norm, 2), minFrN*numel(dataC))); 

data_train = NaN(size(data_norm, 1), minFrN, size(dataC, 1)); % non-NaN pixels x frames per chunk x # train chunks   
data_test = NaN(size(data_norm, 1), minFrN, size(dataC, 1));  % non-NaN pixels x frames per chunk x # test chunks

count_trainsets = 0; 
count_testsets = 0; 
for i = 1:numel(dataC)
    tempDat = data_norm(:, (i-1)*minFrN+1:i*minFrN);
    if i<=size(dataC,1) % because concatenation was columnwise, first column corresponds to train chunks
        count_trainsets = count_trainsets + 1; 
        data_train(:, :, count_trainsets) = tempDat; 
    else
        count_testsets = count_testsets + 1; % second column corresponds to test chunks
        data_test(:, :, count_testsets) = tempDat; 
    end
end
% dff3d = conditionDffMat(data_train(:, :, 1)', nanpxsC{1, 1}, [], [64 64 860]);

%% Save outputs
if ~isempty(save_fn)
    fprintf('\n\tSaving data')
    save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','-v7.3')
    fprintf('\n\tDONE')
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%% HELPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataC_clean = clipDataC_percentile(dataC, lowP, highP)
% clipDataC_percentile
%
% Clips extreme values in a cell array of numeric matrices using global
% percentile thresholds computed across all non-NaN entries in all cells.
%
% INPUTS
%   dataC : cell array of numeric arrays
%   lowP  : lower percentile threshold (default = 0.1)
%   highP : upper percentile threshold (default = 99.9)
%
% OUTPUT
%   dataC_clean : cell array of same size as dataC, with values below/above
%                 the global percentile range clipped to the corresponding
%                 lower/upper threshold
%
% NOTES
%   - Thresholds are computed globally across all cells, not separately
%     within each cell
%   - NaN values are ignored when computing the percentile range
%   - This operation preserves array sizes and avoids introducing new NaNs

if nargin < 2, lowP = 0.1; end
if nargin < 3, highP = 99.9; end

% Collect all non-NaN values across all cells in order to estimate a single
% global clipping range
allVals = [];
for i = 1:numel(dataC)
    if isempty(dataC{i}), continue; end
    v = dataC{i}(:);
    v = v(~isnan(v));
    allVals = [allVals; v]; %#ok<AGROW>
end

% Compute global clipping bounds
lo = prctile(allVals, lowP);
hi = prctile(allVals, highP);

fprintf('Clipping range: [%.3f, %.3f]\n', lo, hi);

% Apply clipping cell-by-cell while preserving original array shapes
dataC_clean = dataC;
for i = 1:numel(dataC)
    if isempty(dataC{i}), continue; end
    dat = dataC{i};
    dat(dat < lo) = lo;
    dat(dat > hi) = hi;
    dataC_clean{i} = dat;
end

end

function data_clean = clipData_percentile(data, lowP, highP)
% clipData_percentile
%
% Clip a numeric data matrix/array based on global percentile thresholds.
%
% Usage:
%   data_clean = clipData_percentile(data)
%   data_clean = clipData_percentile(data, 0.1, 99.9)
%
% Inputs:
%   data  : numeric matrix/array
%   lowP  : lower percentile, default = 0.1
%   highP : upper percentile, default = 99.9
%
% Output:
%   data_clean : data clipped to [lowP, highP] percentile range

if nargin < 2 || isempty(lowP)
    lowP = 0.1;
end

if nargin < 3 || isempty(highP)
    highP = 99.9;
end

if ~isnumeric(data)
    error('Input data must be numeric.');
end

if lowP < 0 || highP > 100 || lowP >= highP
    error('Percentiles must satisfy 0 <= lowP < highP <= 100.');
end

% Collect valid values
allVals = data(:);
allVals = allVals(isfinite(allVals));

if isempty(allVals)
    warning('No finite values found. Returning input unchanged.');
    data_clean = data;
    return
end

% Compute clipping bounds
lo = prctile(allVals, lowP);
hi = prctile(allVals, highP);

fprintf('Clipping range: [%.3f, %.3f]\n', lo, hi);

% Apply clipping while preserving NaNs/Infs as-is unless finite and out of range
data_clean = data;

finiteMask = isfinite(data_clean);

data_clean(finiteMask & data_clean < lo) = lo;
data_clean(finiteMask & data_clean > hi) = hi;

end

function [r_final] = bypassLucric(y, gamma, smt, p_num)
% bypassLucric
%
% Control version of lucric.m that retains the same non-negativity
% enforcement and basic input checks, but intentionally skips the
% Lucy-Richardson deconvolution step.
%
% INPUTS
%   y      : TxN matrix, typically one fluorescence trace per column
%   gamma  : unused here; retained only for interface compatibility with
%            lucric.m
%   smt    : unused here except for reproducing the original input-length
%            validity check
%   p_num  : used only for reproducing the original kernel-length check
%
% OUTPUT
%   r_final : output matrix of same size as y after non-negativity
%             enforcement, with no deconvolution applied
%
% PURPOSE
%   This function is useful as a control for diagnosing whether changes
%   observed after lucric.m arise from:
%     (1) the non-negativity shift itself
%     (2) the deconvolution step
%
% NOTES
%   - MATLAB interprets min(y) columnwise for a 2D matrix, so the line
%     y = y - min(y) shifts each column independently so that its minimum
%     becomes zero
%   - This function is therefore not a no-op; it changes the signal by
%     removing negative values through per-column baseline shifting

% the algorithm works only on strictly non negative input
y = y - min(y);

T = size(y,1);
if T < (p_num*2+2)
    error('Lucy-Richardon requires at least 2*p_num+2 measurement points than its kernel')
end
if smt > T*2-3
    error('Not enough smoothing points are available, choose a smaller smt or a longer trace')
end

% Bypass deconvolution and simply return the shifted signal
r_final = y;

end