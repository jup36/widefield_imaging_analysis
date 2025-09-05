function out = tdr_cross_session_decoder(trainRez, test_tbytDat_hAligned, test_y, varargin)
%TDR_CROSS_SESSION_DECODER  Cross-session generalization of a time-resolved decoder.
%
% out = TDR_CROSS_SESSION_DECODER(trainRez, test_tbytDat_hAligned, test_y, ...)
% applies a time-resolved logistic decoder trained in one session to another
% session (same animal). Uses original-units weights W and intercept b from
% the training session, so no re-standardization is needed.
%
% INPUTS:
%   trainRez               struct from tdr_time_resolved_logistic (.W, .b, .winCtrs, .params[, .Wz])
%   test_tbytDat_hAligned  2xN cell array per trial: {H (K x Tn), t (1 x Tn)} for the TEST session
%   test_y                 Nx1 labels for TEST session (0 = No-Go/CR, 1 = Go/Hit)
%
% OPTIONAL NAME-VALUE:
%   'Epoch'        [t0 t1]  -- override trainRez.params.Epoch for test binning
%   'Win'          scalar   -- override bin width
%   'Step'         scalar   -- override step
%   'minPerClass'  scalar   -- minimum # trials per class in TEST bin to score (default: 5)
%   'DoPlots'      logical  -- plot heatmaps (default: true)
%   'printPlots'   logical  -- save figures as PDF (default: false)
%   'figSaveDir'   char     -- where to save PDFs if printPlots=true
%   'visibleFigs'  logical  -- show or hide figures on screen (default: true)
%   'testRez'      struct   -- optional test session rez (from tdr_time_resolved_logistic) to compute cosine
%   'alignRefEpoch' [a b]   -- reference epoch (sec) for sign alignment when comparing weights (default: [2 4])
%
% OUTPUT struct:
%   out.crossAcc      [nW_train x nW_test] accuracy (train-time x test-time)
%   out.cosineSim     [nW_train x nW_test] cosine similarity of Wz (if testRez provided), else []
%   out.winCtrsTr     1 x nW_train centers used by the TRAIN model
%   out.winCtrsTe     1 x nW_test  centers used to bin the TEST data
%   out.validMask     [N_test x nW_test] which TEST trials contributed per bin
%   out.params        parameters used for test binning/eval
%
% NOTES:
%   - Assumes feature dimensionality K matches between training and test sessions.
%   - Decoding uses linear score s = X_test * W(:,t) + b(t); decision threshold 0 → class 1.
%   - No k-fold CV here (cross-session generalization by design).
%
% Junchol Park, Buschman Lab – 2025


% ---------- Parse inputs ----------
p = inputParser;
p.addParameter('Epoch', [], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('Win',   [], @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Step',  [], @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('minPerClass', 5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('DoPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('printPlots', false, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('figSaveDir', '', @(v)ischar(v)||isstring(v));
p.addParameter('visibleFigs', true, @(v)islogical(v)||ismember(v,[0 1]));
% Across-session cosine options
p.addParameter('testRez', [], @(s)isstruct(s) || isempty(s));
p.addParameter('alignRefEpoch', [2 4], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('alignSigns', false, @(v)islogical(v)||ismember(v,[0 1]));     % NEW
% Naming for saved figures
p.addParameter('trainHeader','', @(v)ischar(v)||isstring(v));                 % NEW
p.addParameter('testHeader', '', @(v)ischar(v)||isstring(v));                 % NEW
p.parse(varargin{:});
prm = p.Results;

% ---------- Pull train info ----------
if ~isfield(trainRez,'W') || ~isfield(trainRez,'b') || ~isfield(trainRez,'winCtrs')
    error('trainRez must contain fields .W, .b, and .winCtrs.');
end
Wtr   = trainRez.W;
btr   = trainRez.b(:)';
ctrTr = trainRez.winCtrs;
[Ktr, nWtr] = size(Wtr);

% ---------- Sanity: feature dimension ----------
Kte = size(test_tbytDat_hAligned{1,1},1);
if Ktr ~= Kte
    error('Feature dimension mismatch: train K=%d, test K=%d.', Ktr, Kte);
end

% ---------- Derive test binning ----------
if ~isfield(trainRez,'params')
    error('trainRez.params is missing; need Epoch/Win/Step to bin the test session.');
end
prm0  = trainRez.params;
Epoch = ternary(~isempty(prm.Epoch), prm.Epoch, prm0.Epoch);
Win   = ternary(~isempty(prm.Win),   prm.Win,   prm0.Win);
Step  = ternary(~isempty(prm.Step),  prm.Step,  prm0.Step);

% ---------- Build test windows ----------
tStart = Epoch(1); tEnd = Epoch(2);
ctrTe  = (tStart + Win/2) : Step : (tEnd - Win/2);
nWte   = numel(ctrTe);
winBoundsTe = [ctrTe(:)-Win/2, ctrTe(:)+Win/2];

% ---------- Assemble test features per bin ----------
N = size(test_tbytDat_hAligned,2);
XwTe = cell(1,nWte); validwTe = cell(1,nWte);
for j = 1:nWte
    lb = winBoundsTe(j,1); ub = winBoundsTe(j,2);
    X = nan(N, Ktr);
    for n = 1:N
        t = test_tbytDat_hAligned{2,n};
        idx = t >= lb & t < ub;
        if any(idx)
            H = test_tbytDat_hAligned{1,n};
            X(n,:) = mean(H(:,idx), 2, 'omitnan')';
        end
    end
    XwTe{j}    = X;
    validwTe{j}= ~any(isnan(X),2);
end

% ---------- Cross-session cross-temporal accuracy ----------
crossAcc  = nan(nWtr, nWte);
validMask = false(N, nWte);
for j = 1:nWte
    validMask(:,j) = validwTe{j};
end

y = test_y(:);
for i = 1:nWtr
    wi = Wtr(:,i);  bi = btr(i);
    if any(isnan(wi)) || isnan(bi), continue; end
    for j = 1:nWte
        valid = validwTe{j};
        if ~any(valid), continue; end
        v1 = sum(valid & (y==1)); v0 = sum(valid & (y==0));
        if v1 < prm.minPerClass || v0 < prm.minPerClass, continue; end
        Xtest  = XwTe{j};
        scores = Xtest(valid,:) * wi + bi;
        yhat   = scores > 0;
        crossAcc(i,j) = mean(yhat == y(valid));
    end
end

% ---------- Across-session cosine similarity (optional) ----------
cosineSim = [];
if ~isempty(prm.testRez) && isfield(trainRez,'Wz') && isfield(prm.testRez,'Wz')
    WzTr = trainRez.Wz; WzTe = prm.testRez.Wz;
    if size(WzTr,1) ~= size(WzTe,1)
        warning('Cannot compute cosineSim: feature dimension mismatch (train vs test Wz).');
    else
        if prm.alignSigns
            WzTrA = align_decoder_signs(WzTr, trainRez.winCtrs, prm.alignRefEpoch);
            WzTeA = align_decoder_signs(WzTe, prm.testRez.winCtrs, prm.alignRefEpoch);
        else
            WzTrA = WzTr; WzTeA = WzTe;
        end
        cosineSim = nan(size(WzTrA,2), size(WzTeA,2));
        for i = 1:size(WzTrA,2)
            wi = WzTrA(:,i); if all(isnan(wi)), continue; end
            for j = 1:size(WzTeA,2)
                wj = WzTeA(:,j); if all(isnan(wj)), continue; end
                mask = isfinite(wi) & isfinite(wj);
                if ~any(mask), continue; end
                ni = norm(wi(mask)); nj = norm(wj(mask));
                if ni==0 || nj==0, continue; end
                cosineSim(i,j) = dot(wi(mask), wj(mask)) / (ni*nj);
            end
        end
    end
end

% ---------- Pack output ----------
out = struct( ...
    'crossAcc',   crossAcc, ...
    'cosineSim',  cosineSim, ...
    'winCtrsTr',  ctrTr, ...
    'winCtrsTe',  ctrTe, ...
    'validMask',  validMask, ...
    'params',     struct('Epoch',Epoch,'Win',Win,'Step',Step,'minPerClass',prm.minPerClass), ...
    'trainHeader', string(prm.trainHeader), ...
    'testHeader',  string(prm.testHeader) );

% ---------- Plotting ----------
if prm.DoPlots
    figVis = ternary(prm.visibleFigs,'on','off');

    % Accuracy heatmap
    hAcc = figure('Name','Cross-session cross-temporal decoding','Visible',figVis);
    imagesc(ctrTr, ctrTe, crossAcc'); axis xy equal tight;
    xlabel('Train time (s)', 'Interpreter','none'); ylabel('Test time (s)', 'Interpreter','none');
    ttl = 'Cross-session cross-temporal decoding accuracy';
    if ~isempty(prm.trainHeader) && ~isempty(prm.testHeader)
        ttl = sprintf('%s\n%s → %s', ttl, prm.trainHeader, prm.testHeader);
    end
    title(ttl, 'Interpreter','none'); colorbar; pbaspect([1 1 1]); clim([0.3 0.9]);
    xline(0,'w:','LineWidth',1.5); xline(2,'w:','LineWidth',1.5);
    yline(0,'w:','LineWidth',1.5); yline(2,'w:','LineWidth',1.5);

    % Cosine heatmap (if computed)
    if ~isempty(cosineSim)
        hCos = figure('Name','Across-session cosine similarity','Visible',figVis);
        imagesc(ctrTr, ctrTe, cosineSim'); axis xy equal tight;
        xlabel('Train time (s)', 'Interpreter','none'); ylabel('Test time (s)', 'Interpreter','none');
        ttl2 = 'Cosine(train W_z, test W_z)';
        if ~isempty(prm.trainHeader) && ~isempty(prm.testHeader)
            ttl2 = sprintf('%s\n%s → %s', ttl2, prm.trainHeader, prm.testHeader);
        end
        title(ttl2, 'Interpreter','none'); colorbar; pbaspect([1 1 1]); clim([-0.3 1]);
        xline(0,'w:','LineWidth',1.5); xline(2,'w:','LineWidth',1.5);
        yline(0,'w:','LineWidth',1.5); yline(2,'w:','LineWidth',1.5);
    else
        hCos = [];
    end

    % Save (yymmdd) & close if invisible
    if prm.printPlots && ~isempty(prm.figSaveDir)
        if ~exist(prm.figSaveDir,'dir'); mkdir(prm.figSaveDir); end
        stamp = char(datetime('now','Format','yyMMdd'));
        base = sprintf('crossSession_%s_%s_%s', string(prm.trainHeader), string(prm.testHeader), stamp);
        print(hAcc, fullfile(prm.figSaveDir, [base '_acc']), '-dpdf','-painters','-bestfit');
        if ~isempty(hCos)
            print(hCos, fullfile(prm.figSaveDir, [base '_cosine']), '-dpdf','-painters','-bestfit');
        end
    end
    if ~prm.visibleFigs
        close(hAcc);
        if ~isempty(hCos), close(hCos); end
    end
end
end

% ===================== Local helpers =====================
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end

function Wz_aligned = align_decoder_signs(Wz, ctr, refEpoch)
Wz_aligned = Wz;
if isempty(Wz) || isempty(ctr) || isempty(refEpoch), return; end
refIdx = find(ctr >= refEpoch(1) & ctr <= refEpoch(2));
if isempty(refIdx), return; end
refVec = mean(Wz(:,refIdx), 2, 'omitnan');
maskR  = isfinite(refVec); nr = norm(refVec(maskR));
if nr==0, return; end
refUnit = zeros(size(refVec)); refUnit(maskR) = refVec(maskR)/nr;
for t = 1:size(Wz,2)
    w = Wz(:,t); mask = isfinite(w) & maskR; if ~any(mask), continue; end
    nw = norm(w(mask)); if nw==0, continue; end
    cosSim = dot(w(mask), refUnit(mask)) / nw;
    if cosSim < 0, Wz_aligned(:,t) = -w; end
end
end