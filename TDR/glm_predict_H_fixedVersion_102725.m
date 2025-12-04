
%% 0) Grab files
filePath = compatiblepath('Z:\Rodent Data\dualImaging_parkj\m1045_jRGECO_GRABda\m1045_122424\task'); 
header = extract_date_animalID_header(filePath); 
keyword_red = '_refitChunks_red_dff_combined.mat';
keyword_beh = '_alignedPupilOrofacial.mat'; 

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath})); 
filePath_H = cell2mat(GrabFiles_sort_trials([header '*' keyword_red], 0, {filePath_matfiles})); 
filePath_B = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles})); 

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

assert(length(tbytDat_hAligned)==length(tbytDat)) % sanity check
N = length(tbytDat); % the number of trials

%% 1) Stack trials to align H
Hs = stack_trials_H(tbytDat_hAligned);  % implement as in your decoder'c', 0.8
Ybig = cell2mat(Hs.Yw'); 
isGoodY = all(isfinite(Ybig), 2);   % drop rows with any NaN/Inf in predictors

%% 2) Build basis for discrete events with pre- and post-event lags
% from stack_trials_H
nW = numel(Hs.winCtrs);
step = Hs.params.Step;   % typically 0.05 s

[basisToneOn.B, basisToneOn.lagsSec] = make_rcos_basis_ortho(9, [0 2], step, ...
    'nonlin','linear', ...
    'c', 0.8, ...   % ~10% of window
    'centerZero', true, ...
    'orthonorm', false);

[basisToneOff.B, basisToneOff.lagsSec] = make_rcos_basis_ortho(9, [0 2], step, ...
    'nonlin','linear', ...
    'c', 0.8, ...   % ~10% of window
    'centerZero', true, ...
    'orthonorm', false);

[basisOutcome.B, basisOutcome.lagsSec] = make_rcos_basis_ortho(9, [-1 1], step, ...
    'nonlin','linear', ...
    'c', 0.8, ...   % ~10% of window
    'centerZero', true, ...
    'orthonorm', false);

[basisLick.B, basisLick.lagsSec] = make_rcos_basis_ortho(9, [-1 1], step, ...
    'nonlin','linear', ...
    'c', 0.8, ...
    'centerZero', true, ...
    'orthonorm', false);

% Optional visualization of basis 
figure('Color','w');
plot(basisToneOff.lagsSec, basisToneOff.B, 'LineWidth',1.2);
    xlabel('Lag (s)'); ylabel('Basis amplitude'); grid on
    title("basis_ToneOff", 'Interpreter','none');
plotBasisKernels(basisToneOn) % visualize with imagesc

%% 3) Build the design matrix (categorical, task, and behavioral variables)
trI = trialTypeInfoAuditoryGngTbytDat(tbytDat); 

%%%%%%%%% discrete variables %%%%%%%%%
% tone onset (Go trials)
toneOn_Go_tbyb = events_to_design(tbytDat, Hs.winBounds, 'evtOn', 'absTime', true, 'trialLogic', trI.goI); % trial x nW
[toneOn_Go_rc, toneOn_Go_rc_names] = convolve_events_with_basis(toneOn_Go_tbyb, basisToneOn, 'toneOnGo'); 

% tone offset (Go trials)
toneOff_Go_tbyb = events_to_design(tbytDat, Hs.winBounds, 'evtOff', 'absTime', true, 'trialLogic', trI.goI); % trial x nW
[toneOff_Go_rc, toneOff_Go_rc_names] = convolve_events_with_basis(toneOff_Go_tbyb, basisToneOff, 'toneOffGo'); 

% tone onset (NoGo trials)
toneOn_NoGo_tbyb = events_to_design(tbytDat, Hs.winBounds, 'evtOn', 'absTime', true, 'trialLogic', trI.nogoI); % trial x nW
[toneOn_NoGo_rc, toneOn_NoGo_rc_names] = convolve_events_with_basis(toneOn_NoGo_tbyb, basisToneOn, 'toneOnNoGo'); 

% tone offset (NoGo trials)
toneOff_NoGo_tbyb = events_to_design(tbytDat, Hs.winBounds, 'evtOff', 'absTime', true, 'trialLogic', trI.nogoI); % trial x nW
[toneOff_NoGo_rc, toneOff_NoGo_rc_names] = convolve_events_with_basis(toneOff_NoGo_tbyb, basisToneOff, 'toneOffNoGo'); 

% licks
periToneLicksVid_tbyb = events_to_design(tbytDat, Hs.winBounds, 'periToneLicksVid', 'absTime', false); % trial x nW
postToneLicksVid_tbyb = events_to_design(tbytDat, Hs.winBounds, 'postToneLicksVid', 'absTime', false); % trial x nW
combinedLicksVid_tbyb = periToneLicksVid_tbyb + postToneLicksVid_tbyb; % trial x nW
[lick_rc, lick_rc_names] = convolve_events_with_basis(combinedLicksVid_tbyb, basisLick, 'lick'); 

% outcome
water_tbyb = events_to_design(tbytDat, Hs.winBounds, 'water', 'absTime', true); % trial x nW
[water_rc, water_rc_names] = convolve_events_with_basis(water_tbyb, basisOutcome, 'water'); 

airpuff_tbyb = events_to_design(tbytDat, Hs.winBounds, 'airpuff', 'absTime', true); % trial x nW
[airpuff_rc, airpuff_rc_names] = convolve_events_with_basis(airpuff_tbyb, basisOutcome, 'airpuff'); 

%%%%%%%%% continuous variables %%%%%%%%%
origCTS = -1:0.02:6; 

% pupil diameter
pupilDiaC = cellfun(@(a) a(1,:), {tbytDat.trPupilDia}, 'UniformOutput', 0); 

X_pupil = timeseries_to_design(pupilDiaC, origCTS, ...
    'Epoch', [-0.9 5], 'Win', 0.1, 'Step', 0.05, ...
    'Method','mean', 'Fill', NaN, 'zscore', true, ...
    'clipExtremes', true, 'zThres', 20);
X_pupil_flat = X_pupil(:); 

% nose tip movement
noseTipC = cellfun(@(a) a(1,:), {tbytDat.trNoseTipEnv}, 'UniformOutput', 0); 

X_nosetip = timeseries_to_design(noseTipC, origCTS, ...
    'Epoch', [-0.9 5], 'Win', 0.1, 'Step', 0.05, ...
    'Method','mean', 'Fill', NaN, 'zscore', true, ...
    'clipExtremes', true, 'zThres', 10);
X_nosetip_flat = X_nosetip(:); 

% whisker movement
whiskerC = cellfun(@(a) a(1,:), {tbytDat.trWhiskerEnv}, 'UniformOutput', 0); 

X_whisker = timeseries_to_design(whiskerC, origCTS, ...
    'Epoch', [-0.9 5], 'Win', 0.1, 'Step', 0.05, ...
    'Method','mean', 'Fill', NaN, 'zscore', true, ...
    'clipExtremes', true, 'zThres', 10);
X_whisker_flat = X_whisker(:); 

% locomotion velocity
locomVelC = cellfun(@(a) a(1,:), {tbytDat.trSpdLocom}, 'UniformOutput', 0); 

X_locomVel = timeseries_to_design(locomVelC, origCTS, ...
    'Epoch', [-0.9 5], 'Win', 0.1, 'Step', 0.05, ...
    'Method','mean', 'Fill', NaN, 'zscore', true, ...
    'clipExtremes', true, 'zThres', 10);
X_locomVel_flat = X_locomVel(:); 

%%%%%%%%% Concatenate design %%%%%%%%%
X_names  = [toneOn_Go_rc_names, toneOn_NoGo_rc_names, toneOff_Go_rc_names, toneOff_NoGo_rc_names, ...
            lick_rc_names, water_rc_names, airpuff_rc_names, ...
            {'pupil'}, {'nosetip'}, {'whisker'}, {'locomVel'}];

X_design = [toneOn_Go_rc, toneOn_NoGo_rc, toneOff_Go_rc, toneOff_NoGo_rc, ...
            lick_rc, water_rc, airpuff_rc, ...
            X_pupil_flat, X_nosetip_flat, X_whisker_flat, X_locomVel_flat];

isGoodX = all(isfinite(X_design), 2);   % drop rows with any NaN/Inf in predictors

X_designVal = X_design(isGoodX & isGoodY, :); 
YbigVal = Ybig(isGoodX & isGoodY, :); 

assert(size(X_designVal, 1)==size(YbigVal, 1)); 

% Standardize X_design
muX = mean(X_designVal, 1, 'omitnan');
sdX = std(X_designVal, 0, 1, 'omitnan');
sdX(~isfinite(sdX) | sdX < 1e-3) = 1e-3;
Xz = (X_designVal - muX) ./ sdX;

%%%%%%%%% Optional diagnosis of colinearity %%%%%%%%%
[fh, R, P, order] = plotDesignCorr(X_design(:,1:9), 'Varnames', X_names(1:9)); 

%% 4) Model fitting
numMotifs = size(YbigVal, 2);
lambdaGrid = logspace(-3, 3, 15);

for k = 1:numMotifs
    yk = YbigVal(:, k);
    keep = all(isfinite(Xz), 2) & isfinite(yk);
    [lambdaBest(k), loss] = pick_lambda_ridge(Xz(keep,:), yk(keep), lambdaGrid, 5);
    beta(:,k) = (Xz(keep,:)'*Xz(keep,:) + lambdaBest(k)*eye(size(Xz,2))) \ (Xz(keep,:)'*yk(keep));
end

%% 5) Explained Variance Analysis with variance partitioning
% Suppose you know the (trial) index of each row (bin-major):
rowTrialIdx = repelem((1:N)', nW);   % N trials, nW bins each (bin-major)
rowTrialIdx = rowTrialIdx(isGoodX & isGoodY); 

% A) Per-motif λ via inner CV (most flexible)
[R2_per, R2_all, det] = cv_global_r2_ridge(Xz, YbigVal, rowTrialIdx, ...
    'outerFolds', 5, ...
    'lambda', [], ...                          % trigger inner CV
    'lambdaGrid', logspace(-3,3,15), ...
    'lambdaMode','perTarget', ...
    'innerFolds', 5, ...
    'standardize', true, ...
    'rng', 1);

% Optional visualization of R2 per motif
plot_R2_per_motif(R2_per)
motifR2Id = find(R2_per > 0.02); 
YbigValR2 = YbigVal(:, motifR2Id); 
motifR2IdC = cell(1, numel(motifR2Id)); 
for kk = 1:sum(motifR2I)
    motifR2IdC{kk} = sprintf("m%d", motifR2Id(kk)); 
end
% A) Recalculate per-motif λ via inner CV (most flexible)
[R2_per_sel, R2_all_sel, det_sel] = cv_global_r2_ridge(Xz, YbigValR2, rowTrialIdx, ...
    'outerFolds', 5, ...
    'lambda', [], ...                          % trigger inner CV
    'lambdaGrid', logspace(-3,3,15), ...
    'lambdaMode','perTarget', ...
    'innerFolds', 5, ...
    'standardize', true, ...
    'rng', 1);

groups = { ...
  struct('name','GoToneOn',   'cols', find(contains(X_names,'toneOnGo'))), ...
  struct('name','NoGoToneOn', 'cols', find(contains(X_names,'toneOnNoGo'))), ...
  struct('name','ToneOffGo',  'cols', find(contains(X_names,'toneOffGo'))), ...
  struct('name','ToneOffNoGo','cols', find(contains(X_names,'toneOffNoGo'))), ...
  struct('name','Lick',       'cols', find(contains(X_names,'lick_rc'))), ...
  struct('name','Water',      'cols', find(contains(X_names,'water_rc'))), ...
  struct('name','Airpuff',    'cols', find(contains(X_names,'airpuff_rc'))), ...
  struct('name','Pupil',      'cols', find(strcmp(X_names,'pupil'))), ...
  struct('name','Nose',       'cols', find(strcmp(X_names,'nosetip'))), ...
  struct('name','Whisker',    'cols', find(strcmp(X_names,'whisker'))), ...
  struct('name','Locom',      'cols', find(strcmp(X_names,'locomVel'))) ...
};

glmEvRez = variance_partition_timeShuffle( ...
          Xz, YbigVal, lambdaBest, ...
          'KFold', 5, ...
          'RowIndex', rowTrialIdx, ...
          'Groups', groups, ...
          'ShuffleMode', 'within_trial', ...
          'RngSeed', 1, ...
          'Verbose', true);

% Examples:
% glmEvRez.cvR2_full           -> K x 1
% glmEvRez.cvR2_unique         -> P x K
% glmEvRez.cvR2_unique_grp     -> G x K (group-wise ΔR²)
% glmEvRez.cvR2_upper_grp      -> G x K

%% 6) Visualize beta weights (optional)
plot_beta_with_labels(beta(1:36, :), X_names(1:36), 22)

[bw.S, bw.order, bw.hDendro, bw.hSim] = beta_cosine_map(beta); 

figSaveDir = fullfile(filePath, "Figure"); 
if exist(figSaveDir, "dir")~=7
    mkdir(figSaveDir)
end

timestamp = char(datetime('now', 'Format', 'MMddyy'));  % format: mmddyy_hhmm
print(bw.hDendro, fullfile(figSaveDir, ['glm_beta_dendrogram', '_', header, '_' timestamp]), '-dpdf', '-painters', '-bestfit')
print(bw.hSim, fullfile(figSaveDir, ['glm_beta_cosine_similarity', '_', header, '_' timestamp]), '-dpdf', '-painters', '-bestfit')

%% 7) Save Model
% -------------- pack output --------------
glmRez = struct();
glmRez.beta            = beta;          % P x K
glmRez.X_names         = X_names;       % 1 x P 
glmRez.lambdaBest      = lambdaBest;    % 1 x K
glmRez.X_design        = X_design;      % (NxnW) x P 
glmRez.Ybig            = Ybig;          % (NxnW) x K 
glmRez.isGoodX         = isGoodX;       % (NxnW) x 1
glmRez.isGoodY         = isGoodY;       % (NxnW) x 1 

saveDir = fullfile(filePath_matfiles, [header, '_glmRez_', timestamp, '.mat']);  
save(fullfile(saveDir), 'glmRez', 'glmEvRez', '-v7.3')










