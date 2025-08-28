function out = tdr_cross_time_decoder_analysis(figSaveDir, tbytDat_hAligned, y, varargin)
%TDR_CROSS_TIME_DECODER_ANALYSIS  Analyze generalization of decoding across time.
%
% This function trains time-resolved logistic regression classifiers using
% TDR (targeted dimensionality reduction) to decode binary trial labels
% (e.g., hit vs. CR) from neural population activity aligned to behavioral events.
% 
% For each time window, a classifier is trained and tested on all time bins to 
% produce:
%   - Cosine similarity matrix of decoder weights (Wz) across time
%   - Cross-temporal decoding accuracy matrix
%
% Plots and saves PDF figures with cue epoch (0–2s) highlighted.
%
% INPUTS:
%   figSaveDir        - directory to save output figures (PDFs)
%   tbytDat_hAligned  - 1x2 cell array per trial: {H matrix, time vector}
%   y                 - Nx1 binary trial outcome labels (0 or 1)
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'Epoch'           - [start, end] of analysis window (default: [-0.5 4])
%   'Win'             - width of time window (default: 0.1 sec)
%   'Step'            - step size between windows (default: 0.05 sec)
%   'Lambda'          - ridge penalty for logistic regression (default: 1)
%   'minSD'           - minimum feature std for z-scoring (default: 1e-3)
%   'minPerClass'     - min trials per class to train (default: 5)
%   'DoPlots'         - flag to display plots (default: true)
%   'printPlots'      - flag to save plots as PDFs (default: true)
%
% OUTPUT (struct):
%   out.Wz            - raw decoder weights (features x time bins)
%   out.Wz_aligned    - sign-aligned decoder weights (if enabled)
%   out.b             - bias terms per time bin
%   out.b_aligned     - sign-aligned bias terms (if enabled)
%   out.cosineSim     - cosine similarity matrix of decoder weights
%   out.crossDecodingAcc - decoding accuracy when training/testing across time
%   out.winCtrs       - center times of decoding windows
%
% DEPENDENCIES:
%   align_decoder_signs (included as a nested helper function)
%
% Junchol Park, Buschman Lab – 2025


% ---- Parameters ----
p = inputParser;
p.addParameter('Epoch', [-0.5 4], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('Win', 0.100, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Step', 0.050, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Lambda', 1, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('minSD', 1e-3, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('minPerClass', 5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('DoPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('printPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('visibleFigs', true, @(v)islogical(v)||ismember(v,[0 1]));
p.parse(varargin{:});
prm = p.Results;

% ---- Time window setup ----
tStart = prm.Epoch(1); tEnd = prm.Epoch(2);
win = prm.Win; step = prm.Step;
ctr = (tStart + win/2) : step : (tEnd - win/2);
nW = numel(ctr);
winBounds = [ctr(:)-win/2, ctr(:)+win/2];

% ---- Feature assembly ----
N = size(tbytDat_hAligned, 2);
K = size(tbytDat_hAligned{1,1},1);
Xw = cell(1,nW); validw = cell(1,nW);

for i = 1:nW
    lb = winBounds(i,1); ub = winBounds(i,2);
    X = nan(N, K);
    for n = 1:N
        t = tbytDat_hAligned{2,n};
        idx = t >= lb & t < ub;
        if any(idx)
            H = tbytDat_hAligned{1,n};
            X(n,:) = mean(H(:,idx), 2, 'omitnan')';
        end
    end
    Xw{i} = X;
    validw{i} = ~any(isnan(X),2);
end

% ---- Global z-scoring stats ----
Xstack = cat(1, Xw{:});
muG = mean(Xstack, 1, 'omitnan');
sdG = std(Xstack, 0, 1, 'omitnan');
sdG(sdG < prm.minSD) = prm.minSD;

% ---- Train classifiers at each time ----
Wz = nan(K, nW); b = nan(1,nW);
validTrials = false(N, nW);
for i = 1:nW
    X = Xw{i};
    valid = validw{i};
    v_hit = valid & (y(:)==1);
    v_cr  = valid & (y(:)==0);
    if sum(v_hit) < prm.minPerClass || sum(v_cr) < prm.minPerClass
        continue
    end
    Xz = (X(valid,:) - muG) ./ sdG;
    yv = y(valid);
    mdl = fitclinear(Xz, yv, 'Learner','logistic', 'Regularization','ridge', ...
        'Lambda', prm.Lambda, 'Solver','lbfgs', 'ClassNames',[0,1]);
    Wz(:,i) = mdl.Beta;
    b(i) = mdl.Bias;
    validTrials(valid,i) = true;
end

% Align decoder weights based on response period (e.g. 2–4s)
%[Wz_aligned, flipSign] = align_decoder_signs(Wz, ctr, [2 4]);
%b_aligned = b;
%b_aligned(flipSign) = -b_aligned(flipSign); % sign flip must be applied here too
Wz_aligned = Wz; b_aligned = b;

% ---- Cosine similarity and decoding accuracy across time ----
cosineSim = nan(nW, nW);
decAcc    = nan(nW, nW);

for i = 1:nW
    wi = Wz_aligned(:,i);
    if all(isnan(wi)), continue; end

    for j = 1:nW
        wj = Wz_aligned(:,j);
        if all(isnan(wj)), continue; end

        % --- robust cosine similarity (ignore NaNs; guard zero norms)
        mask = isfinite(wi) & isfinite(wj);
        if ~any(mask), continue; end
        ni = norm(wi(mask)); nj = norm(wj(mask));
        if ni==0 || nj==0, continue; end
        cosineSim(i,j) = dot(wi(mask), wj(mask)) / (ni * nj);

        % --- cross-time decoding: train at i (wi,b_aligned(i)), test at j
        Xtest = Xw{j};
        valid = validw{j} & validTrials(:,i);
        if ~any(valid), continue; end

        Xz_test = (Xtest(valid,:) - muG) ./ sdG;
        y_test  = y(valid);

        scores = Xz_test * wi + b_aligned(i);   % use aligned bias

        yhat   = scores > 0;
        decAcc(i,j) = mean(yhat == y_test);
    end
    fprintf("Completed cross-time evaluation for win#%d/%d\n", i, nW)
end

% ---- Pack output ----
out = struct('Wz', Wz, 'Wz_aligned', Wz_aligned, 'b', b, 'b_aligned', b_aligned, ...
            'cosineSim', cosineSim, 'crossDecodingAcc', decAcc, 'winCtrs', ctr);

% ---- Plots ----
if prm.DoPlots
    figVis = 'on';
    if ~prm.visibleFigs, figVis = 'off'; end

    h_cos = figure('Visible', figVis);
    imagesc(ctr, ctr, cosineSim); axis xy equal tight;
    title('Cosine Similarity of Decoder Weights'); xlabel('Train Time'); ylabel('Test Time'); colorbar;
    pbaspect([1 1 1]); clim([-0.3 1])
    % Add white dotted cue lines (0–2 s)
    xline(0,  'w:','LineWidth',1.5);
    xline(2,  'w:','LineWidth',1.5);
    yline(0,  'w:','LineWidth',1.5);
    yline(2,  'w:','LineWidth',1.5);

    h_dcd = figure('Visible', figVis); 
    imagesc(ctr, ctr, decAcc); axis xy equal tight;
    title('Cross-temporal Decoding Accuracy'); xlabel('Train Time'); ylabel('Test Time'); colorbar;
    pbaspect([1 1 1]); clim([0.3 0.9])
    xline(0,  'w:','LineWidth',1.5);
    xline(2,  'w:','LineWidth',1.5);
    yline(0,  'w:','LineWidth',1.5);
    yline(2,  'w:','LineWidth',1.5);
end

% Print (Optional)
if prm.printPlots
    if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
    header = extract_date_animalID_header(figSaveDir);
    timestampStr = char(datetime('now','Format','MMddyy_HHmmss'));
    figSaveName_cos = sprintf('across_time_w_cosSimilarity_%s_%s', header, timestampStr);
    print(h_cos, fullfile(figSaveDir, figSaveName_cos),'-dpdf','-painters','-bestfit')
    figSaveName_dcd = sprintf('across_time_decodeAccuracy_%s_%s', header, timestampStr);
    print(h_dcd, fullfile(figSaveDir, figSaveName_dcd),'-dpdf','-painters','-bestfit')
end

if ~prm.visibleFigs
    close(h_cos); close(h_dcd);
end

%% % ---- Helper function ----
function [Wz_aligned, flipSign] = align_decoder_signs(Wz, ctr, refEpoch)
% Aligns decoder weights Wz across time bins based on a reference epoch
%   Wz       - K x T matrix of decoder weights (K: features, T: time bins)
%   ctr      - 1 x T vector of window centers (in seconds)
%   refEpoch - [start, end] reference period for sign alignment (e.g., [2, 4])

    T = size(Wz,2);
    Wz_aligned = Wz;
    flipSign = false(1, T);

    % Identify reference bins
    refIdx = find(ctr >= refEpoch(1) & ctr <= refEpoch(2));
    if isempty(refIdx)
        warning('No reference bins found in specified refEpoch.');
        return
    end

    % Compute reference vector (mean over ref bins)
    refVec = mean(Wz(:,refIdx), 2, 'omitnan');
    refVecNorm = norm(refVec(~isnan(refVec)));
    if refVecNorm == 0
        warning('Reference vector has zero norm; skipping alignment.');
        return
    end
    refUnit = refVec / refVecNorm;

    % Align each time bin
    for tt = 1:T
        w = Wz(:,tt);
        if all(isnan(w)), continue; end
        val = isfinite(w) & isfinite(refUnit);
        if ~any(val), continue; end
        cosSim = dot(w(val), refUnit(val)) / (norm(w(val)) * norm(refUnit(val)));
        if cosSim < 0
            Wz_aligned(:,tt) = -w;
            flipSign(tt) = true;
        end
    end
end

end
