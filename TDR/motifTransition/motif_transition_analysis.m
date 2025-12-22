% filePath = 'Z:\Rodent Data\dualImaging_parkj\m1045_jRGECO_GRABda\m1045_122424\task'; 
% fileKeyword = '_refitChunksCurated_red_dff_combined.mat';

%% 0) Grab files
header      = extract_date_animalID_header(filePath);
keyword_beh = '_alignedPupilOrofacial.mat';

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath}));
filePath_H        = cell2mat(GrabFiles_sort_trials([header '*' fileKeyword], 0, {filePath_matfiles}));
filePath_B        = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

assert(length(tbytDat_hAligned)==length(tbytDat)) % sanity check
N = length(tbytDat); % # trials

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);

%% 1) Stack trials to align H
Hs   = stack_trials_H(tbytDat_hAligned, 'zscore', true);
nBin500ms = round(0.5/Hs.params.Step); 

[XcorrMat_sess, lags] = computeMotifXcorr(Hs.Y3, nBin500ms); 

% --- choose lag range for "future" ---
posLagMask = (lags > 0 & lags <= nBin500ms);   % example: first 5 positive bins

% --- directional transition-like matrix ---
A_dir = squeeze(mean(XcorrMat_sess(:,:,posLagMask), 3));   % K x K

% figure('Color','w');
% imagesc(A_dir);
% clim([-0.1 0.1])
% axis xy; colorbar;



