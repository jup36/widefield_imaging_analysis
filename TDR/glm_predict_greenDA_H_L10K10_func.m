function glm_predict_greenDA_H_L10K10_func(filePath, filePath_H, filePath_B)
%% GLM over green DA motif time courses with robust design construction
%
% SYNOPSIS
%   glm_predict_greenDA_H_L10K10_func(filePath, filePath_H, filePath_B)
%
% DESCRIPTION
%   Fits ridge-regularized GLMs to predict per-motif green DA temporal
%   activity (H) from task events and behavioral covariates in a single
%   session.
%
% PIPELINE
%   1) Load session data:
%        - combined green DA refit file containing tbytDat_hAligned
%        - behavioral tbytDat file
%   2) Stack trials into bin-major matrices.
%   3) Build event regressors using raised-cosine basis functions.
%   4) Add continuous behavioral regressors only if present.
%   5) Standardize X and z-score Y.
%   6) Fit ridge GLM per motif.
%   7) Compute trial-wise CV R².
%   8) Run variance partitioning.
%   9) Save GLM output and figures with greenDA-specific names.
%
% INPUTS
%   filePath : full path to session task folder
%              e.g. ...\m1045_122424\task
%
%   filePath_H : full path to combined green DA refit output
%                must contain tbytDat_hAligned
%
%   filePath_B : full path to behavioral tbytDat file
%                must contain tbytDat
%
% OUTPUT
%   Saves:
%       <Matfiles>/<header>_glmRez_greenDA_L10K10_<timestamp>.mat
%
% Author: Junchol Park
% Updated: green DA version

%% 0) Validate inputs and grab files

filePath_matfiles = fullfile(filePath, 'Matfiles');

if exist(filePath_matfiles, 'dir') ~= 7
    error("Matfiles folder doesn't exist:\n%s", filePath_matfiles);
end

if iscell(filePath_B)
    filePath_B = filePath_B{1};
end

if exist(filePath_B, 'file') ~= 2
    error("Behavior data file doesn't exist:\n%s", filePath_B);
end

if iscell(filePath_H)
    filePath_H = filePath_H{1};
end

if exist(filePath_H, 'file') ~= 2
    error("Combined green DA refit file doesn't exist:\n%s", filePath_H);
end

if ~contains(filePath_H, 'greenDA') && ~contains(filePath_H, 'green_dff')
    warning('filePath_H does not appear to be a greenDA file:\n%s', filePath_H);
end

header = extract_date_animalID_header(filePath_B);

fprintf('\n============================================================\n');
fprintf('GLM prediction for green DA H motifs\n');
fprintf('============================================================\n');
fprintf('Session/task path: %s\n', filePath);
fprintf('H/refit file:      %s\n', filePath_H);
fprintf('Behavior file:     %s\n', filePath_B);
fprintf('Header:            %s\n', header);

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

assert(length(tbytDat_hAligned) == length(tbytDat), ...
    'tbytDat_hAligned and tbytDat have different lengths.');

N = length(tbytDat);

emptyHtrials = cellfun(@isempty, tbytDat_hAligned(1, :));

fprintf('tbytDat_hAligned empty trials: %d/%d\n', ...
    sum(emptyHtrials), numel(emptyHtrials));

if any(emptyHtrials)
    fprintf('Empty trial indices:\n');
    disp(find(emptyHtrials));
end

%% 1) Stack trials to align H

Hs   = stack_trials_H(tbytDat_hAligned, 'zscore', false);
Ybig = cell2mat(Hs.Yw');         % [(N*nW) x K], raw motif activity

isGoodY = all(isfinite(Ybig), 2); %#ok<NASGU>

%% 2) Basis for discrete events

nW   = numel(Hs.winCtrs);
step = Hs.params.Step;

[basisToneOn.B,  basisToneOn.lagsSec]  = make_rcos_basis_ortho(9, [0 2], step, ...
    'nonlin', 'linear', 'c', 0.8, 'centerZero', true, 'orthonorm', false);

[basisToneOff.B, basisToneOff.lagsSec] = make_rcos_basis_ortho(9, [0 2], step, ...
    'nonlin', 'linear', 'c', 0.8, 'centerZero', true, 'orthonorm', false);

[basisOutcome.B, basisOutcome.lagsSec] = make_rcos_basis_ortho(9, [-1 1], step, ...
    'nonlin', 'linear', 'c', 0.8, 'centerZero', true, 'orthonorm', false);

[basisLick.B, basisLick.lagsSec] = make_rcos_basis_ortho(9, [-1 1], step, ...
    'nonlin', 'linear', 'c', 0.8, 'centerZero', true, 'orthonorm', false);

%% 3) Design matrix: discrete event regressors

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);

% Tone on/off, Go/NoGo
toneOn_Go_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'evtOn', 'absTime', true, 'trialLogic', trI.goI);

[toneOn_Go_rc, toneOn_Go_rc_names] = convolve_events_with_basis( ...
    toneOn_Go_tbyb, basisToneOn, 'toneOnGo');

toneOff_Go_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'evtOff', 'absTime', true, 'trialLogic', trI.goI);

[toneOff_Go_rc, toneOff_Go_rc_names] = convolve_events_with_basis( ...
    toneOff_Go_tbyb, basisToneOff, 'toneOffGo');

toneOn_NoGo_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'evtOn', 'absTime', true, 'trialLogic', trI.nogoI);

[toneOn_NoGo_rc, toneOn_NoGo_rc_names] = convolve_events_with_basis( ...
    toneOn_NoGo_tbyb, basisToneOn, 'toneOnNoGo');

toneOff_NoGo_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'evtOff', 'absTime', true, 'trialLogic', trI.nogoI);

[toneOff_NoGo_rc, toneOff_NoGo_rc_names] = convolve_events_with_basis( ...
    toneOff_NoGo_tbyb, basisToneOff, 'toneOffNoGo');

% Licks
periToneLicksVid_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'periToneLicksVid', 'absTime', false);

postToneLicksVid_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'postToneLicksVid', 'absTime', false);

combinedLicksVid_tbyb = periToneLicksVid_tbyb + postToneLicksVid_tbyb;

[lick_rc, lick_rc_names] = convolve_events_with_basis( ...
    combinedLicksVid_tbyb, basisLick, 'lick');

% Outcomes
water_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'water', 'absTime', true);

[water_rc, water_rc_names] = convolve_events_with_basis( ...
    water_tbyb, basisOutcome, 'water');

airpuff_tbyb = events_to_design( ...
    tbytDat, Hs.winBounds, 'airpuff', 'absTime', true);

[airpuff_rc, airpuff_rc_names] = convolve_events_with_basis( ...
    airpuff_tbyb, basisOutcome, 'airpuff');

X_blocks = {toneOn_Go_rc, toneOn_NoGo_rc, ...
            toneOff_Go_rc, toneOff_NoGo_rc, ...
            lick_rc, water_rc, airpuff_rc};

X_names = [toneOn_Go_rc_names, toneOn_NoGo_rc_names, ...
           toneOff_Go_rc_names, toneOff_NoGo_rc_names, ...
           lick_rc_names, water_rc_names, airpuff_rc_names];

%% 3b) Continuous variables, conditionally included

origCTS = -1:0.02:6;

[X_blocks, X_names] = try_add_cont( ...
    tbytDat, 'trPupilDia', 'pupil', 20, origCTS, X_blocks, X_names);

[X_blocks, X_names] = try_add_cont( ...
    tbytDat, 'trNoseTipEnv', 'nosetip', 10, origCTS, X_blocks, X_names);

[X_blocks, X_names] = try_add_cont( ...
    tbytDat, 'trWhiskerEnv', 'whisker', 10, origCTS, X_blocks, X_names);

[X_blocks, X_names] = try_add_cont( ...
    tbytDat, 'trSpdLocom', 'locomVel', 10, origCTS, X_blocks, X_names);

%% 3c) Concatenate design

X_design = cat(2, X_blocks{:});

X_designVal = X_design;
YbigVal = Ybig;

%% 3d) Z-score Y

Ymu = mean(YbigVal, 1, 'omitnan');
Ysd = std(YbigVal, 0, 1, 'omitnan');
Ysd(~isfinite(Ysd) | Ysd < 1e-6) = 1;

Yz = (YbigVal - Ymu) ./ Ysd;
Yz(~isfinite(Yz)) = 0;

%% 3e) Robust standardization and NaN handling for X

nPred = size(X_designVal, 2);
Xz = nan(size(X_designVal));
muX = nan(1, nPred);
sdX = nan(1, nPred);

for j = 1:nPred
    xj = X_designVal(:, j);

    mu = mean(xj, 'omitnan');
    sd = std(xj, 0, 'omitnan');

    if ~isfinite(sd) || sd < 1e-3
        sd = 1e-3;
    end

    zj = (xj - mu) ./ sd;
    zj(~isfinite(zj)) = 0;

    Xz(:, j) = zj;
    muX(j) = mu;
    sdX(j) = sd;
end

% Drop ultra-sparse predictors based on finite fraction before z-scoring.
finiteFrac = mean(isfinite(X_designVal), 1);
dropCols = finiteFrac < 0.25;

if any(dropCols)
    fprintf('Dropping %d sparse predictor columns.\n', sum(dropCols));

    Xz(:, dropCols) = [];
    X_names(dropCols) = [];
    muX(dropCols) = [];
    sdX(dropCols) = [];
end

%% 4) Ridge GLM fitting per motif

numMotifs = size(Yz, 2);
lambdaGrid = logspace(-3, 3, 15);

beta = nan(size(Xz, 2), numMotifs);
lambdaBest = nan(1, numMotifs);

for k = 1:numMotifs
    yk_raw = Yz(:, k);
    keep = isfinite(yk_raw);

    [lambdaBest(k), ~] = pick_lambda_ridge( ...
        Xz(keep, :), yk_raw(keep), lambdaGrid, 5);

    XtX = Xz(keep, :)' * Xz(keep, :);
    Xty = Xz(keep, :)' * yk_raw(keep);

    beta(:, k) = (XtX + lambdaBest(k) * eye(size(Xz, 2))) \ Xty;
end

%% 5) Descriptive stats for estimated and observed H

Hest_statC = motif_event_estimates( ...
    beta, X_names, muX, sdX, Xz, ...
    toneOn_Go_tbyb, toneOff_Go_tbyb, ...
    toneOn_NoGo_tbyb, toneOff_NoGo_tbyb, ...
    basisToneOn, basisToneOff, ...
    Hs.winCtrs, trI, ...
    'ReturnTrials', false);

HsZ = stack_trials_H(tbytDat_hAligned, 'zscore', true);
Horg_statC = descriptiveH(HsZ.Y3, HsZ.winCtrs, trI);

%% 6) Decoder: Go vs NoGo from green DA motif expression

figSaveDirLogit = fullfile(filePath, "Figure", "logitReg_greenDA");

if exist(figSaveDirLogit, "dir") ~= 7
    mkdir(figSaveDirLogit);
end

decBins = tdr_time_resolved_logistic_bins( ...
    figSaveDirLogit, HsZ.Y3, HsZ.winCtrs, trI.goI, ...
    'Lambda', 1, ...
    'minPerClass', 5, ...
    'DoPlots', true, ...
    'printPlots', true, ...
    'visibleFigs', false);

%% 7) CV R² with trial-wise folds

rowTrialIdx = repelem((1:N)', nW);

[R2_per, R2_all, R2_det] = cv_global_r2_ridge( ...
    Xz, Yz, rowTrialIdx, ...
    'outerFolds', 5, ...
    'lambda', [], ...
    'lambdaGrid', logspace(-3, 3, 15), ...
    'lambdaMode', 'perTarget', ...
    'innerFolds', 5, ...
    'standardize', false, ...
    'rng', 1);

%% 8) Variance partitioning

groups = {};

groups{end+1} = struct( ...
    'name', 'GoToneOn', ...
    'cols', find(contains(X_names, 'toneOnGo')));

groups{end+1} = struct( ...
    'name', 'NoGoToneOn', ...
    'cols', find(contains(X_names, 'toneOnNoGo')));

groups{end+1} = struct( ...
    'name', 'ToneOffGo', ...
    'cols', find(contains(X_names, 'toneOffGo')));

groups{end+1} = struct( ...
    'name', 'ToneOffNoGo', ...
    'cols', find(contains(X_names, 'toneOffNoGo')));

groups{end+1} = struct( ...
    'name', 'Lick', ...
    'cols', find(contains(X_names, 'lick_rc')));

groups{end+1} = struct( ...
    'name', 'Water', ...
    'cols', find(contains(X_names, 'water_rc')));

groups{end+1} = struct( ...
    'name', 'Airpuff', ...
    'cols', find(contains(X_names, 'airpuff_rc')));

maybeNames = {'pupil', 'nosetip', 'whisker', 'locomVel'};

for nm = maybeNames
    idx = find(strcmp(X_names, nm{1}));
    if ~isempty(idx)
        groups{end+1} = struct('name', nm{1}, 'cols', idx); %#ok<AGROW>
    end
end

groups = groups(~cellfun(@(g) isempty(g.cols), groups));

glmEvRez = variance_partition_timeShuffle( ...
    Xz, Yz, lambdaBest, rowTrialIdx, ...
    'KFold', 5, ...
    'Groups', groups, ...
    'ShuffleMode', 'within_trial', ...
    'RngSeed', 1, ...
    'Verbose', true);

%% 9) Visualize beta similarity

[bw.S, bw.order, bw.hDendro, bw.hSim] = beta_cosine_map( ...
    beta, 'Title', 'greenDA β cosine (ridge)');

figSaveDir = fullfile(filePath, "Figure");

if exist(figSaveDir, "dir") ~= 7
    mkdir(figSaveDir);
end

timestamp = char(datetime('now', 'Format', 'MMddyy'));

print(bw.hDendro, fullfile(figSaveDir, ...
    ['glm_beta_dendrogram_greenDA_', header, '_', timestamp]), ...
    '-dpdf', '-painters', '-bestfit');

print(bw.hSim, fullfile(figSaveDir, ...
    ['glm_beta_cosine_similarity_greenDA_', header, '_', timestamp]), ...
    '-dpdf', '-painters', '-bestfit');

close([bw.hDendro, bw.hSim]);

%% 10) Save GLM output

glmRez = struct();

glmRez.beta       = beta;
glmRez.X_names    = X_names;
glmRez.lambdaBest = lambdaBest;

glmRez.X_design = X_design;
glmRez.Ybig     = Ybig;
glmRez.Yz       = Yz;
glmRez.Ymu      = Ymu;
glmRez.Ysd      = Ysd;

glmRez.group  = groups;
glmRez.R2_per = R2_per;
glmRez.R2_all = R2_all;
glmRez.det    = R2_det;

glmRez.Hest_statC = Hest_statC;
glmRez.Horg_statC = Horg_statC;

glmRez.muX = muX;
glmRez.sdX = sdX;

glmRez.decBins = decBins;

glmRez.filePath    = filePath;
glmRez.filePath_H  = filePath_H;
glmRez.filePath_B  = filePath_B;
glmRez.header      = header;
glmRez.signalType  = 'greenDA';
glmRez.lagKtag     = 'L10K10';
glmRez.timestamp   = timestamp;
glmRez.emptyHtrials = emptyHtrials;

savePath = fullfile(filePath_matfiles, ...
    [header, '_glmRez_greenDA_L10K10_', timestamp, '.mat']);

save(savePath, 'glmRez', 'glmEvRez', '-v7.3');

fprintf('\nCompleted successful green DA GLM fitting run for %s\n', header);
fprintf('Saved GLM result:\n%s\n', savePath);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_blocks, X_names] = try_add_cont(tbytDat, varField, outName, zThres, origCTS, X_blocks, X_names)
% try_add_cont
%
% Adds one continuous variable to the GLM design matrix only if the field
% exists in tbytDat.

if ~isfield(tbytDat, varField)
    fprintf('Continuous field %s not found. Skipping %s.\n', varField, outName);
    return
end

targetLen = numel(origCTS);

sigC = cell(1, numel(tbytDat));

for ii = 1:numel(tbytDat)

    if isfield(tbytDat, varField) && ~isempty(tbytDat(ii).(varField))
        v = tbytDat(ii).(varField);
        v = v(1, :);
    else
        v = nan(1, targetLen);
    end

    sigC{ii} = v;
end

Xc = timeseries_to_design( ...
    sigC, origCTS, ...
    'Epoch', [-0.9 5], ...
    'Win', 0.1, ...
    'Step', 0.05, ...
    'Method', 'mean', ...
    'Fill', NaN, ...
    'zscore', false, ...
    'clipExtremes', true, ...
    'zThres', zThres);

Xc = Xc(:);

if ~isempty(Xc) && size(Xc, 1) == size(X_blocks{1}, 1)
    X_blocks{end+1} = Xc; %#ok<AGROW>
    X_names{end+1} = outName; %#ok<AGROW>

    fprintf('Added continuous predictor: %s from %s\n', outName, varField);
else
    warning(['Continuous predictor %s from %s was not added because size did not match.\n' ...
             '  size(Xc,1) = %d\n' ...
             '  expected   = %d'], ...
             outName, varField, size(Xc, 1), size(X_blocks{1}, 1));
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HstatC = descriptiveH(HsY3d, Htime, trI) %#ok<INUSD>
% descriptiveH
%
% Computes mean and SEM of motif activity separately for Go and NoGo trials.
%
% INPUTS
%   HsY3d : [N x K x T] motif activity
%   Htime : [1 x T] time vector, retained for interface consistency
%   trI   : trial type struct with goI and nogoI
%
% OUTPUT
%   HstatC : 1 x K cell array of descriptive statistics

HstatC = cell(1, size(HsY3d, 2));

for k = 1:size(HsY3d, 2)

    goMat = squeeze(HsY3d(trI.goI, k, :));
    nogoMat = squeeze(HsY3d(trI.nogoI, k, :));

    [mean_go_H, ~, sem_go_H] = meanstdsem(goMat);
    [mean_nogo_H, ~, sem_nogo_H] = meanstdsem(nogoMat);

    Hstat = struct();

    Hstat.motifIdx = k;

    Hstat.mean = struct( ...
        'periCueGoH', mean_go_H, ...
        'periCueNoGoH', mean_nogo_H);

    Hstat.sem = struct( ...
        'periCueGoH', sem_go_H, ...
        'periCueNoGoH', sem_nogo_H);

    HstatC{1, k} = Hstat;
end

end