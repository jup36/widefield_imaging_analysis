function [data_norm, nanpxs, data_train, data_test] = ProcessAndSplitDataGreenDA(fn,save_fn,parameter_class)
%Camden MacDowell - timeless
%filters, normalizes, and splits data in fn into training and test test
%fn can be the full file path or a stack and opt structure. 
if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

gp = loadobj(feval(parameter_class)); 

%% SUPER UGLY CONTIGENCIES DUE TO LEGACY DATA
if ischar(fn) %load data
    warning('Camden you may want to check how things are transposed')
    fprintf('\n\tLoading data')
    %load data
    temp = load(fn);
    data = temp.dff;
    opts = temp.opts;
    clear temp;
elseif isstruct(fn)
    data = fn.dff;
    opts = fn.opts;
else
    data = fn;
    opts = gp;
end
%%

%condition data and remove nan pxls
[x,y,z] = size(data);   
[data,nanpxs] = conditionDffMat(data); %nanpxs are the same on each iteration so fine to overwrite

% Clip extreme values globally across all chunks to reduce the influence of
% rare artifactual outliers before subsequent processing
data = clipData_percentile(data, 0.2, 99.8); % Mask extreme values

%% Apply bypassed Lucy-Richardson preprocessing pixelwise
% This step preserves the non-negativity enforcement from lucric.m
% but skips the actual deconvolution. Each column is a pixel trace over time.
for px = 1:size(data,2)
    data(:,px) = bypassLucric(data(:,px), gp.d_gamma, gp.d_smooth, gp.d_kernel);
end

%Denoise with PCA (removed banded pixels)
if gp.w_pca_denoise
    data = conditionDffMat(data,nanpxs,[], [x,y,z]);
    data = DenoisePCA(data);
    [data,~] = conditionDffMat(data);
end

%normalize to 0 to 1 
fprintf('\n\tPerforming %s normalization to %d value', gp.w_normalization_method, gp.w_norm_val);
switch gp.w_normalization_method
    case 'pixelwise' %each between x and xth pixel intensity
        data_norm = NaN(size(data));
        for px = 1:size(data,2)
            data_norm(:,px) = (data(:,px))/(prctile(data(:,px),gp.w_norm_val));
        end             
    case 'full' %normalize using the percentile of the maximum         
        data_norm = data/prctile(data(data>eps),gp.w_norm_val);          
    case 'bounded'
        data_norm = (data)/(gp.w_norm_val(2)); %normalize between zero and the upper bound     
    case 'none'
        data_norm = data;
    otherwise
        error('Unknown normalization method. Check general params')
end

%transpose (fpCNMF operates rowwise)
data_norm = data_norm';

%Chunk (since MU is suboptimal for >4500 timepoints. As many as possible
num_chunks = floor(z/(gp.w_chunk_dur*opts.fps));
if ~isEven(num_chunks)% need even number
    num_chunks = num_chunks-1; 
end

%remove the remainder and reshape into chunks
data_trim = data_norm(:,1:end-mod(z,num_chunks*(gp.w_chunk_dur*opts.fps)));

data_chunked = cell(1,num_chunks);
for i = 1:num_chunks 
    data_chunked{i} = data_trim(:,1+(i-1)*(gp.w_chunk_dur*opts.fps):(i*(gp.w_chunk_dur*opts.fps)));
end

%alternate testing and training
data_train = cat(3,data_chunked{1:2:num_chunks});
data_test = cat(3,data_chunked{2:2:num_chunks});

%save off the data in the scratch directory and the nanpxs
if ~isempty(save_fn)
    fprintf('\n\tSaving data')
    save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','-v7.3')
    fprintf('\n\tDONE')
end


%% for saving off figure use; 
% [path, name] = fileparts(save_fn);
% saveCurFigs(gcf,'-dpng',sprintf('registration_%s',fn),[path filesep name, '_ProcessAndSplitFigures'],0); %close all;


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
% Control version of lucric.m that implements the same preprocessing
% steps, but skips the Lucy-Richardson deconvolution.
%
% Inputs:
%   y      - TxN matrix, each column is a fluorescence trace over time
%   gamma  - unused here; included for interface compatibility
%   smt    - unused here; included for interface compatibility
%   p_num  - used only for the same input-length check as lucric
%
% Output:
%   r_final - TxN matrix after enforcing non-negativity, without deconvolution
%
% This function is useful for testing whether changes seen after lucric
% are due to:
%   (1) the non-negativity shift
%   versus
%   (2) the deconvolution itself.

% the algorithm works only on strictly non negative input
y = y - min(y);

T = size(y,1);
if T < (p_num*2+2)
    error('Lucy-Richardon requires at least 2*p_num+2 measurement points than its kernel')
end
if smt > T*2-3
    error('Not enough smoothing points are available, choose a smaller smt or a longer trace')
end

% bypass deconvolution: just return the shifted signal
r_final = y;

end










