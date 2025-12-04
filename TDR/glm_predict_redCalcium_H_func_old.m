function glm_predict_redCalcium_H_func_OLD(filePath, fileKeyword)
%% DEPRECATED (currently using 'glm_predict_redCalcium_H_func.m')
%% GLM over motif timecourses with robust design construction
% glm_predict_redCalcium_H_func(filePath, fileKeyword)
%
% SYNOPSIS
%   Fits ridge-regularized GLMs to predict per-motif temporal activity (H)
%   from task events and behavioral covariates in a single session. The
%   design matrix is built from (i) convolved discrete events (tone on/off,
%   outcomes, licks) and (ii) continuous signals (pupil, nose, whisker,
%   locomotion) that are **included only if present** in `tbytDat`. Missing
%   continuous fields are skipped automatically; names and columns stay
%   consistent throughout the pipeline.
%
% PIPELINE
%   1) Load session data:
%        • _refitChunks_*_dff_combined.mat  → motif activity aligned per trial (H)
%        • _tbytDat_alignedPupilOrofacial.mat → trial structure & behavioral signals
%   2) Stack trials into bin-major matrices (Y) and define analysis window.
%   3) Build event regressors by convolving impulses with raised-cosine bases:
%        toneOn (Go/NoGo), toneOff (Go/NoGo), water, airpuff, licks.
%   4) Conditionally add continuous regressors (if fields exist in tbytDat):
%        pupil (trPupilDia), nose (trNoseTipEnv), whisker (trWhiskerEnv),
%        locomotion (trSpdLocom); each windowed/stepped and z-scored.
%   5) Standardize X; fit ridge GLM per motif with λ chosen by inner CV.
%   6) Compute **trial-wise** outer-CV R² (per motif and global).
%   7) Variance partitioning via temporal shuffling (unique & upper-bound ΔR²)
%      over user-level groups inferred from X_names.
%   8) Save model artifacts (β, λ, R², groups, names, X/Y) and export
%      β-similarity visualizations (dendrogram, cosine heatmap).
%
% INPUTS
%   filePath    : char/str, full path to the session “task” folder
%                 (e.g., '...\m1045_122424\task')
%   fileKeyword : char/str, filename suffix for the calcium file to load
%                 (e.g., '_refitChunks_red_dff_combined.mat' or curated variant)
%
% OUTPUTS (saved to disk)
%   <Matfiles>/<header>_glmRez_<timestamp>.mat containing:
%     glmRez.beta        : [P x K] ridge coefficients (P predictors, K motifs)
%     glmRez.X_names     : 1xP cellstr of predictor names actually used
%     glmRez.lambdaBest  : [1 x K] best λ per motif
%     glmRez.X_design    : [(N*nW) x P] (pre-filter) design matrix
%     glmRez.Ybig        : [(N*nW) x K] motif targets (pre-filter)
%     glmRez.isGoodX/Y   : logical masks used to select valid rows
%     glmRez.group       : grouping struct array used for variance partition
%     glmRez.R2_per      : [1 x K] outer-CV R² per motif (trial-wise folds)
%     glmRez.R2_all      : scalar outer-CV R² pooled across motifs
%     glmRez.det         : details from CV routine (fold R², λ, scalers, …)
%   Also saves `glmEvRez` (variance partition results).
%   Figures (PDF): β dendrogram and cosine-similarity heatmap in task/Figure.
%
% ROBUSTNESS / DESIGN INCLUSION
%   • Continuous predictors are appended only if their fields exist in tbytDat.
%     (No placeholder columns; X_names mirrors actual design.)
%   • Discrete event regressors are always built from trial events.
%   • Rows with any NaN/Inf in X or Y are dropped consistently.
%
% CROSS-VALIDATION
%   • Outer CV splits by trials (not samples) using `rowTrialIdx`.
%   • Inner CV (per motif) selects λ over a log-spaced grid.
%
% NOTES
%   • Hs (motif activity) is taken as produced by `stack_trials_H`.
%   • Basis functions: raised-cosine with specified windows and step sizes.
%   • Standardization: predictors z-scored on the analysis subset.
%
% EXAMPLE
%   fp  = 'Z:\Rodent Data\dualImaging_parkj\m1045_jRGECO_GRABda\m1045_122424\task';
%   key = '_refitChunks_red_dff_combined.mat';
%   glm_predict_redCalcium_H_func(fp, key);


%% 0) Grab files
header = extract_date_animalID_header(filePath);
keyword_beh = '_alignedPupilOrofacial.mat';

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath}));
filePath_H = cell2mat(GrabFiles_sort_trials([header '*' fileKeyword], 0, {filePath_matfiles}));
filePath_B = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

assert(length(tbytDat_hAligned)==length(tbytDat)) % sanity check
N = length(tbytDat); % # trials

%% 1) Stack trials to align H
Hs   = stack_trials_H(tbytDat_hAligned, 'zscore', false);
Ybig = cell2mat(Hs.Yw');
isGoodY = all(isfinite(Ybig), 2);

%% 2) Basis for discrete events
nW   = numel(Hs.winCtrs);
step = Hs.params.Step;   % typically 0.05 s

[basisToneOn.B,  basisToneOn.lagsSec]  = make_rcos_basis_ortho(9, [0 2], step,'nonlin','linear','c',0.8,'centerZero',true,'orthonorm',false);
[basisToneOff.B, basisToneOff.lagsSec] = make_rcos_basis_ortho(9, [0 2], step,'nonlin','linear','c',0.8,'centerZero',true,'orthonorm',false);
[basisOutcome.B, basisOutcome.lagsSec] = make_rcos_basis_ortho(9, [-1 1],step,'nonlin','linear','c',0.8,'centerZero',true,'orthonorm',false);
[basisLick.B,    basisLick.lagsSec]    = make_rcos_basis_ortho(9, [-1 1],step,'nonlin','linear','c',0.8,'centerZero',true,'orthonorm',false);

%% 3) Design matrix (discrete first)
trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);

% tone on/off (Go/NoGo)
toneOn_Go_tbyb    = events_to_design(tbytDat, Hs.winBounds, 'evtOn',  'absTime', true, 'trialLogic', trI.goI);
[toneOn_Go_rc, toneOn_Go_rc_names] = convolve_events_with_basis(toneOn_Go_tbyb, basisToneOn, 'toneOnGo');

toneOff_Go_tbyb   = events_to_design(tbytDat, Hs.winBounds, 'evtOff', 'absTime', true, 'trialLogic', trI.goI);
[toneOff_Go_rc, toneOff_Go_rc_names] = convolve_events_with_basis(toneOff_Go_tbyb, basisToneOff, 'toneOffGo');

toneOn_NoGo_tbyb  = events_to_design(tbytDat, Hs.winBounds, 'evtOn',  'absTime', true, 'trialLogic', trI.nogoI);
[toneOn_NoGo_rc, toneOn_NoGo_rc_names] = convolve_events_with_basis(toneOn_NoGo_tbyb, basisToneOn, 'toneOnNoGo');

toneOff_NoGo_tbyb = events_to_design(tbytDat, Hs.winBounds, 'evtOff', 'absTime', true, 'trialLogic', trI.nogoI);
[toneOff_NoGo_rc, toneOff_NoGo_rc_names] = convolve_events_with_basis(toneOff_NoGo_tbyb, basisToneOff, 'toneOffNoGo');

% licks (peri + post)
periToneLicksVid_tbyb = events_to_design(tbytDat, Hs.winBounds, 'periToneLicksVid', 'absTime', false);
postToneLicksVid_tbyb = events_to_design(tbytDat, Hs.winBounds, 'postToneLicksVid', 'absTime', false);
combinedLicksVid_tbyb = periToneLicksVid_tbyb + postToneLicksVid_tbyb;
[lick_rc, lick_rc_names] = convolve_events_with_basis(combinedLicksVid_tbyb, basisLick, 'lick');

% outcomes
water_tbyb   = events_to_design(tbytDat, Hs.winBounds, 'water',   'absTime', true);
[water_rc,   water_rc_names]   = convolve_events_with_basis(water_tbyb, basisOutcome, 'water');

airpuff_tbyb = events_to_design(tbytDat, Hs.winBounds, 'airpuff', 'absTime', true);
[airpuff_rc, airpuff_rc_names] = convolve_events_with_basis(airpuff_tbyb, basisOutcome, 'airpuff');

% Start with discrete blocks
X_blocks = {toneOn_Go_rc, toneOn_NoGo_rc, toneOff_Go_rc, toneOff_NoGo_rc, ...
            lick_rc, water_rc, airpuff_rc};
X_names  = [toneOn_Go_rc_names, toneOn_NoGo_rc_names, toneOff_Go_rc_names, toneOff_NoGo_rc_names, ...
            lick_rc_names, water_rc_names, airpuff_rc_names];

%% Robust continuous variables (conditionally include if present)
origCTS = -1:0.02:6;

% Helper to try add one continuous var
    function try_add_cont(varField, outName, zThres)
        if isfield(tbytDat, varField)
            % extract 1×T per trial (allow empties)
            sigC = cellfun(@(a) a(1,:), {tbytDat.(varField)}, 'UniformOutput', 0);
            Xc = timeseries_to_design(sigC, origCTS, ...
                    'Epoch', [-0.9 5], 'Win', 0.1, 'Step', 0.05, ...
                    'Method','mean', 'Fill', NaN, 'zscore', true, ...
                    'clipExtremes', true, 'zThres', zThres);
            Xc = Xc(:); % column
            % guard: only append if non-empty and has rows
            if ~isempty(Xc) && size(Xc,1) == size(X_blocks{1},1)
                X_blocks{end+1} = Xc; %#ok<AGROW>
                X_names{end+1}  = outName; %#ok<AGROW>
            end
        end
    end

% Add each continuous predictor only if present
try_add_cont('trPupilDia',   'pupil',    20);
try_add_cont('trNoseTipEnv', 'nosetip',  10);
try_add_cont('trWhiskerEnv', 'whisker',  10);
try_add_cont('trSpdLocom',   'locomVel', 10);

% Concatenate design
X_design = cat(2, X_blocks{:});
isGoodX  = all(isfinite(X_design), 2);

X_designVal = X_design(isGoodX & isGoodY, :);
YbigVal     = Ybig(isGoodX & isGoodY, :);
assert(size(X_designVal,1)==size(YbigVal,1));

% Standardize X
muX = mean(X_designVal, 1, 'omitnan');
sdX = std(X_designVal, 0, 1, 'omitnan');
sdX(~isfinite(sdX) | sdX < 1e-3) = 1e-3;
Xz = (X_designVal - muX) ./ sdX;

%% 4) Model fitting (ridge per motif)
numMotifs   = size(YbigVal, 2);
lambdaGrid  = logspace(-3, 3, 15);
beta        = nan(size(Xz,2), numMotifs);
lambdaBest  = nan(1, numMotifs);

for k = 1:numMotifs
    yk = YbigVal(:, k);
    keep = all(isfinite(Xz), 2) & isfinite(yk);
    [lambdaBest(k), ~] = pick_lambda_ridge(Xz(keep,:), yk(keep), lambdaGrid, 5);
    beta(:,k) = (Xz(keep,:)'*Xz(keep,:) + lambdaBest(k)*eye(size(Xz,2))) \ (Xz(keep,:)'*yk(keep));
end

%% 5) CV R² with trial-wise folds
rowTrialIdx = repelem((1:N)', nW);
rowTrialIdx = rowTrialIdx(isGoodX & isGoodY);

[R2_per, R2_all, R2_det] = cv_global_r2_ridge(Xz, YbigVal, rowTrialIdx, ...
    'outerFolds', 5, ...
    'lambda', [], ...
    'lambdaGrid', logspace(-3,3,15), ...
    'lambdaMode','perTarget', ...
    'innerFolds', 5, ...
    'standardize', true, ...
    'rng', 1);

%% 6) Variance partitioning (groups built from present X_names)
groups = {};
groups{end+1} = struct('name','GoToneOn',   'cols', find(contains(X_names,'toneOnGo')));
groups{end+1} = struct('name','NoGoToneOn', 'cols', find(contains(X_names,'toneOnNoGo')));
groups{end+1} = struct('name','ToneOffGo',  'cols', find(contains(X_names,'toneOffGo')));
groups{end+1} = struct('name','ToneOffNoGo','cols', find(contains(X_names,'toneOffNoGo')));
groups{end+1} = struct('name','Lick',       'cols', find(contains(X_names,'lick_rc')));
groups{end+1} = struct('name','Water',      'cols', find(contains(X_names,'water_rc')));
groups{end+1} = struct('name','Airpuff',    'cols', find(contains(X_names,'airpuff_rc')));
maybeNames = {'pupil','nosetip','whisker','locomVel'};
for nm = maybeNames
    idx = find(strcmp(X_names, nm{1}));
    if ~isempty(idx)
        groups{end+1} = struct('name', nm{1}, 'cols', idx);
    end
end
% remove empty groups (no cols)
groups = groups(~cellfun(@(g) isempty(g.cols), groups));

glmEvRez = variance_partition_timeShuffle( ...
    Xz, YbigVal, lambdaBest, rowTrialIdx, ...
    'KFold', 5, ...
    'Groups', groups, ...
    'ShuffleMode', 'within_trial', ...
    'RngSeed', 1, ...
    'Verbose', true);

%% 7) Visualize β (invisible figs)
[bw.S, bw.order, bw.hDendro, bw.hSim] = beta_cosine_map(beta, 'Title','β cosine (ridge)');
figSaveDir = fullfile(filePath, "Figure");
if exist(figSaveDir, "dir")~=7, mkdir(figSaveDir); end
timestamp = char(datetime('now', 'Format', 'MMddyy'));
print(bw.hDendro, fullfile(figSaveDir, ['glm_beta_dendrogram_', header, '_', timestamp]), '-dpdf', '-painters', '-bestfit');
print(bw.hSim,    fullfile(figSaveDir, ['glm_beta_cosine_similarity_', header, '_', timestamp]), '-dpdf', '-painters', '-bestfit');
close([bw.hDendro, bw.hSim]);

%% 8) Save
glmRez = struct();
glmRez.beta       = beta;
glmRez.X_names    = X_names;
glmRez.lambdaBest = lambdaBest;
glmRez.X_design   = X_design;
glmRez.Ybig       = Ybig;
glmRez.isGoodX    = isGoodX;
glmRez.isGoodY    = isGoodY;
glmRez.group      = groups;
glmRez.R2_per     = R2_per;
glmRez.R2_all     = R2_all;
glmRez.det        = R2_det;

saveDir = fullfile(filePath_matfiles, [header, '_glmRez_', timestamp, '.mat']);
save(fullfile(saveDir), 'glmRez', 'glmEvRez', '-v7.3');
end