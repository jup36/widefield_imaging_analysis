function RefitCuratedBasisMotifs_JP(fn, basis_dir, chunk, parameter_class, save_dir, varargin)
% RefitCuratedBasisMotifs_JP
%
% Refit curated spatiotemporal basis motifs to one data chunk.
%
% This version uses fpCNMF_refit, which keeps W fixed and only fits H.
%
% Inputs:
%   fn              : processed file containing data_train, data_test, nanpxs
%   basis_dir       : .mat file containing W_basis
%   chunk           : chunk index
%   parameter_class : e.g. 'general_params_dual'
%   save_dir        : output directory
%
% Optional:
%   'dateStr'             : string used in output filename
%   'H_init_method'       : 'projection', 'constant', or 'provided'
%   'H_init_prctile_cap'  : percentile cap for projection initializer

%% Optional inputs
p = inputParser;
addParameter(p, 'dateStr', datestr(now, 'mmddyy'), @(x) ischar(x) || isstring(x));
addParameter(p, 'H_init_method', 'projection', @(x) ischar(x) || isstring(x));
addParameter(p, 'H_init_prctile_cap', 99.5, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});

dateStr = char(p.Results.dateStr);
H_init_method = char(p.Results.H_init_method);
H_init_prctile_cap = p.Results.H_init_prctile_cap;

%% Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));
end

if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

gp = loadobj(feval(parameter_class));

%% Load basis motifs
fprintf('\n\tLoading basis motifs\n');

basisS = load(basis_dir, 'W_basis');
W_full = basisS.W_basis;

%% Load train/test chunk data
fprintf('\n\tLoading data\n');

dataS = load(fn, 'data_test', 'data_train', 'nanpxs');

if chunk < 1 || chunk > size(dataS.data_test, 3)
    error('Invalid chunk index %d. data_test has %d chunks.', chunk, size(dataS.data_test, 3));
end

data_train = dataS.data_train(:,:,chunk);
data_test  = dataS.data_test(:,:,chunk);
nanpxs = dataS.nanpxs;

[~, name] = fileparts(fn);
name = [name, sprintf('_refitChunkCurated_%s_%d', dateStr, chunk)];

%% Recondition, smooth, and flatten
W = W_full;

if numel(gp.smt_kernel) > 2

    fprintf('\n\tAutofitting smoothing value\n');

    gp.smt_kernel = AutoFitSmoothingLevel(cat(2, data_train, data_test), nanpxs, W, gp);

    data_train = reshape( ...
        SpatialGaussian(conditionDffMat(data_train', nanpxs), gp.smt_kernel), ...
        gp.pixel_dim(1)*gp.pixel_dim(2), size(data_train, 2));

    data_test = reshape( ...
        SpatialGaussian(conditionDffMat(data_test', nanpxs), gp.smt_kernel), ...
        gp.pixel_dim(1)*gp.pixel_dim(2), size(data_test, 2));

elseif numel(gp.smt_kernel) == 2

    fprintf('\n\tUsing set smoothing value\n');

    data_train = reshape( ...
        SpatialGaussian(conditionDffMat(data_train', nanpxs), gp.smt_kernel), ...
        gp.pixel_dim(1)*gp.pixel_dim(2), size(data_train, 2));

    data_test = reshape( ...
        SpatialGaussian(conditionDffMat(data_test', nanpxs), gp.smt_kernel), ...
        gp.pixel_dim(1)*gp.pixel_dim(2), size(data_test, 2));

else

    fprintf('\n\tNo smoothing\n');

    data_train = reshape( ...
        conditionDffMat(data_train', nanpxs), ...
        gp.pixel_dim(1)*gp.pixel_dim(2), size(data_train, 2));

    data_test = reshape( ...
        conditionDffMat(data_test', nanpxs), ...
        gp.pixel_dim(1)*gp.pixel_dim(2), size(data_test, 2));
end

%% Only work on shared pixels
fprintf('\n\tRestricting to shared pixels\n');

[~, bad_pxl] = SharedPixels(W, cat(2, data_test, data_train));

data_train(bad_pxl,:) = [];
data_test(bad_pxl,:)  = [];
W(bad_pxl,:,:)        = [];

W_fixed_masked = W;

fprintf('\n\tMasked W size: %s\n', mat2str(size(W_fixed_masked)));
fprintf('\tTrain data size: %s\n', mat2str(size(data_train)));
fprintf('\tTest data size: %s\n', mat2str(size(data_test)));

%% ------------------------------------------------------------------------
% Train refit
% -------------------------------------------------------------------------

fprintf('\n\tRefitting train chunk %d\n', chunk);

[w, H, stats_internal, refitInfo] = fpCNMF_refit(data_train, ...
    'non_penalized_iter', gp.non_penalized_iter, ...
    'penalized_iter', gp.penalized_iter_refit, ...
    'speed', 'fast', ...
    'verbose', 0, ...
    'lambda', 0, ...
    'ortho_H', gp.ortho_H, ...
    'W', W_fixed_masked, ...
    'sparse_H', 0, ...
    'H_init_method', H_init_method, ...
    'H_init_prctile_cap', H_init_prctile_cap);

stats_refit = CNMF_Stats(w, H, data_train, 0);
stats_refit.smoothingkernel = gp.smt_kernel;
stats_refit.refitInfo = refitInfo;
stats_refit.stats_internal = stats_internal;

residuals = data_train - tensor_convolve(w, H);

save(fullfile(save_dir, [name 'train.mat']), ...
    'w', 'H', 'stats_refit', 'bad_pxl', 'residuals', ...
    'W_fixed_masked', 'dateStr', 'chunk', 'H_init_method', 'H_init_prctile_cap', ...
    '-v7.3');

fprintf('\n\tSaved train refit: %s\n', fullfile(save_dir, [name 'train.mat']));

%% ------------------------------------------------------------------------
% Test refit
% -------------------------------------------------------------------------

fprintf('\n\tRefitting test chunk %d\n', chunk);

[w, H, stats_internal, refitInfo] = fpCNMF_refit(data_test, ...
    'non_penalized_iter', gp.non_penalized_iter, ...
    'penalized_iter', gp.penalized_iter_refit, ...
    'speed', 'fast', ...
    'verbose', 0, ...
    'lambda', 0, ...
    'ortho_H', gp.ortho_H, ...
    'W', W_fixed_masked, ...
    'sparse_H', 0, ...
    'H_init_method', H_init_method, ...
    'H_init_prctile_cap', H_init_prctile_cap);

stats_refit = CNMF_Stats(w, H, data_test, 0);
stats_refit.smoothingkernel = gp.smt_kernel;
stats_refit.refitInfo = refitInfo;
stats_refit.stats_internal = stats_internal;

residuals = data_test - tensor_convolve(w, H);

save(fullfile(save_dir, [name 'test.mat']), ...
    'w', 'H', 'stats_refit', 'bad_pxl', 'residuals', ...
    'W_fixed_masked', 'dateStr', 'chunk', 'H_init_method', 'H_init_prctile_cap', ...
    '-v7.3');

fprintf('\n\tSaved test refit: %s\n', fullfile(save_dir, [name 'test.mat']));

fprintf('\n\tEnd of a successful curated fixed-W refit run\n');

end