function [xcorrPosLagMat, xcorrPosLagMat_go, xcorrPosLagMat_nogo, ...
    XcorrMat_sess, XcorrMat_sess_go, XcorrMat_sess_nogo] = ...
    motifH_crosscorr_posLag_func(filePath, fileKeyword, varargin)
% MOTIFH_CROSSCORR_POSLAG_FUNC
%   Compute motif–motif cross-correlogram summary (mean positive-lag xcorr)
%   for one session, optionally allowing the user to specify the post-lag
%   duration in seconds (default = 0.5 s).
%
% INPUTS
%   filePath       : full path to session "task" folder
%   fileKeyword    : suffix for H-file (e.g. '_refitChunksCurated_red_dff_combined.mat')
%
% NAME–VALUE PAIRS
%   'postLagPeriod' : scalar, in seconds (default: 0.5)
%
% OUTPUTS
%   xcorrPosLagMat      : [K x K] summary for ALL trials
%   xcorrPosLagMat_go   : [K x K] summary for Go trials
%   xcorrPosLagMat_nogo : [K x K] summary for No-Go trials
%
% NOTES
%   Positive lags (future correlations) summarize how motif j tends to
%   follow motif i within the postLagPeriod window.

%% -------------------- Parse inputs --------------------
p = inputParser;
p.addParameter('postLagPeriod', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
postLagPeriod_sec = p.Results.postLagPeriod;

%% -------------------- 0) Grab files --------------------
header      = extract_date_animalID_header(filePath);
keyword_beh = '_alignedPupilOrofacial.mat';

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath}));
filePath_H        = cell2mat(GrabFiles_sort_trials([header '*' fileKeyword], 0, {filePath_matfiles}));
filePath_B        = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);
assert(size(tbytDat_hAligned, 2) == numel(tbytDat), 'Mismatch in trial counts.');

%% -------------------- 1) Stack trials --------------------
Hs = stack_trials_H(tbytDat_hAligned, 'zscore', true);
[~, K, ~] = size(Hs.Y3);

stepSec = Hs.params.Step;
nBin = round(postLagPeriod_sec / stepSec);   % convert seconds → bins

%% -------------------- 2) Compute cross-correlograms --------------------
% ALL TRIALS
[XcorrMat_sess, lags] = computeMotifXcorr(Hs.Y3, nBin);

% GO TRIALS
if any(trI.goI)
    X_go = Hs.Y3(trI.goI, :, :);
    XcorrMat_sess_go = computeMotifXcorr(X_go, nBin);
else
    XcorrMat_sess_go = nan(K, K, 2*nBin + 1);
end

% NOGO TRIALS
if any(trI.nogoI)
    X_ng = Hs.Y3(trI.nogoI, :, :);
    XcorrMat_sess_nogo = computeMotifXcorr(X_ng, nBin);
else
    XcorrMat_sess_nogo = nan(K, K, 2*nBin + 1);
end

%% -------------------- 3) Extract positive-lag portion --------------------
posLagMask = (lags > 0 & lags <= nBin);

xcorrPosLagMat      = squeeze(mean(XcorrMat_sess(:, :,      posLagMask), 3));
xcorrPosLagMat_go   = squeeze(mean(XcorrMat_sess_go(:, :,   posLagMask), 3));
xcorrPosLagMat_nogo = squeeze(mean(XcorrMat_sess_nogo(:, :, posLagMask), 3));

end


%% ------------------------------------------------------------------------
%% Helper: full cross-correlogram
function [XcorrMat, lags] = computeMotifXcorr(HsY3, maxLag)

if nargin < 2 || isempty(maxLag), maxLag = 10; end
if maxLag < 0 || maxLag ~= round(maxLag)
    error('maxLag must be a nonnegative integer.');
end

[N, K, T] = size(HsY3);
L    = 2*maxLag + 1;
lags = -maxLag:maxLag;

XcorrMat = zeros(K, K, L);

for i = 1:K
    Xi = reshape(HsY3(:, i, :), [], 1);
    Xi(~isfinite(Xi)) = 0;

    for j = 1:K
        Xj = reshape(HsY3(:, j, :), [], 1);
        Xj(~isfinite(Xj)) = 0;

        XcorrMat(i, j, :) = xcorr(Xi, Xj, maxLag, 'coeff');
    end
end
end
