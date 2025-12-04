function out = tdr_time_resolved_logistic_bins(figSaveDir, HsY3d, Htime, y, varargin)
%TDR_TIME_RESOLVED_LOGISTIC_BINS Time-resolved TDR via logistic regression (per time bin).
%
% out = tdr_time_resolved_logistic_bins(figSaveDir, HsY3d, Htime, y, ...)
%
% INPUTS
%   figSaveDir : directory to save output figures (string/char)
%
%   HsY3d      : [N x K x T] motif activity (e.g., HsZ.Y3)
%                N = #trials, K = #motifs, T = #time bins
%
%   Htime      : [1 x T] time vector (e.g., Hs.winCtrs)
%
%   y          : [N x 1] labels (0 = NoGo, 1 = Go)
%
% Name–Value pairs
%   'Lambda'       : ridge regularization strength (default 1)
%   'minSD'        : min SD for variance flooring in unscaling (default 1e-3)
%   'minPerClass'  : min trials per class to fit model (default 5)
%   'DoPlots'      : whether to generate plots (default true)
%   'printPlots'   : whether to save plots as PDF (default true)
%   'cdWeightCscale' : color scale for weight heatmap (default [-0.15 0.15])
%   'dPrimeYlim'   : y-limits for d' plot (default [0.2 2])
%   'projCdYlim'   : y-limits for projection plot (default [-1.5 2])
%   'visibleFigs'  : show figures on screen (default true)
%
% OUTPUT
%   out : struct with fields
%       .Wz       : [K x T] standardized weights (per bin)
%       .W        : [K x T] unscaled weights (per bin)
%       .b        : [1 x T] intercepts (per bin, in original-units space)
%       .proj     : [N x T] projections per trial & time bin
%       .dprime   : [1 x T] d' across time
%       .motifRank: [T x K] motif indices sorted by |weight| per bin
%       .params   : struct of parameters
%       .cvAcc    : [1 x T] 5-fold CV decoding accuracy per bin
%       .time     : [1 x T] copy of Htime
%
% NOTES
%   • This is *bin-wise* decoding: at each time bin t, X = H(:, :, t) (N x K).
%   • Logistic regression (ridge) is fit separately for each bin.
%   • Sign alignment: Go (y==1) trials are forced to have more positive
%     projections than NoGo on average at each bin.
%   • Global μ/σ over all bins are used for standardization (like your
%     windowed TDR code), so weights are comparable across time.
%

% ---------------- parameters ----------------
p = inputParser;
p.addParameter('Lambda', 1, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('minSD', 1e-3, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('minPerClass', 5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('DoPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('printPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('cdWeightCscale', [-0.15 0.15], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('dPrimeYlim', [0.2 2], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('projCdYlim', [-1.5 2], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('visibleFigs', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('rng', 1, @(v)isscalar(v));
p.parse(varargin{:});
prm = p.Results;

% ---------------- basics ----------------
[N, K, T] = size(HsY3d);
assert(numel(Htime)==T, 'Htime must be length T = size(HsY3d,3).');
assert(numel(y)==N, 'Label vector y must have length N = size(HsY3d,1).');

y = y(:);                       % ensure column
goMask   = (y == 1);
nogoMask = (y == 0);

% ---------------- global stats across all trials & time bins ----------------
% Stack H into (N*T) x K for global μ/σ
Hstack = reshape(HsY3d, N*T, K);    % (N*T) x K
muG = mean(Hstack, 1, 'omitnan');   % 1 x K
sdG = std( Hstack, 0, 1, 'omitnan');
sdG(~isfinite(sdG) | sdG < prm.minSD) = prm.minSD;

% ---------------- alloc output ----------------
Wz   = nan(K, T);     % standardized weights
W    = nan(K, T);     % unscaled weights
b    = nan(1, T);     % intercepts
proj = nan(N, T);     % projections
dprime = nan(1, T);   % d'
cvAcc  = nan(1, T);   % decoding accuracy
motifRank = nan(T, K);

rng(prm.rng);

for t = 1:T
    % features at this time bin: N x K
    X = squeeze(HsY3d(:, :, t));   % [N x K]

    % valid trials: finite features
    valid = all(isfinite(X), 2);
    if ~any(valid), continue; end

    y_valid = y(valid);
    X_valid = X(valid, :);

    % check class counts
    maskGo   = (y_valid == 1);
    maskNoGo = (y_valid == 0);
    if sum(maskGo) < prm.minPerClass || sum(maskNoGo) < prm.minPerClass
        continue;
    end

    % global z-scoring (same μ/σ across all bins)
    Xz_valid = (X_valid - muG) ./ sdG;
    Xz_valid(~isfinite(Xz_valid)) = 0;

    % ------------- 5-fold CV decoding accuracy -------------
    cv = cvpartition(y_valid, 'KFold', 5);
    acc_fold = nan(cv.NumTestSets, 1);

    for kfold = 1:cv.NumTestSets
        trIdx = training(cv, kfold);
        teIdx = test(cv, kfold);

        mdl_k = fitclinear(Xz_valid(trIdx,:), y_valid(trIdx), ...
            'Learner','logistic', 'Regularization','ridge', ...
            'Lambda', prm.Lambda, 'Solver','lbfgs', ...
            'ClassNames',[0,1]);

        yhat_k = predict(mdl_k, Xz_valid(teIdx,:));
        acc_fold(kfold) = mean(yhat_k == y_valid(teIdx));
    end

    cvAcc(t) = mean(acc_fold, 'omitnan');

    % ------------- final model using all valid trials -------------
    mdl = fitclinear(Xz_valid, y_valid, ...
        'Learner','logistic', 'Regularization','ridge', ...
        'Lambda', prm.Lambda, 'Solver','lbfgs', ...
        'ClassNames',[0,1]);

    w_z = mdl.Beta;      % K x 1 (standardized)
    b_t = mdl.Bias;

    % store standardized weights
    Wz(:,t) = w_z;

    % unscale to original units (same as in your windowed code)
    w_unscaled = w_z ./ sdG(:);            % K x 1
    b_unz      = b_t - muG * (w_z ./ sdG(:));

    % projections for all trials (keep NaNs for invalid ones)
    s = nan(N,1);
    s(valid) = X_valid * w_unscaled + b_unz;

    % sign alignment: Go > NoGo
    go_s   = s(goMask);
    nogo_s = s(nogoMask);
    if mean(go_s, 'omitnan') < mean(nogo_s, 'omitnan')
        w_unscaled = -w_unscaled;
        w_z        = -w_z;
        s          = -s;
        b_unz      = -b_unz;
    end

    W(:,t)    = w_unscaled;
    proj(:,t) = s;
    b(t)      = b_unz;

    % d' at this bin
    go_s   = s(goMask);
    nogo_s = s(nogoMask);
    m1 = mean(go_s, 'omitnan');
    m0 = mean(nogo_s, 'omitnan');
    v1 = nanvar(go_s);
    v0 = nanvar(nogo_s);
    dprime(t) = (m1 - m0) / sqrt(0.5 * (v1 + v0) + eps);

    % motif rank by |weight|
    [~, ord] = sort(abs(w_unscaled),'descend');
    motifRank(t,:) = ord(:)';
end

% ---------------- pack output ----------------
out = struct();
out.Wz      = Wz;
out.W       = W;
out.b       = b;
out.proj    = proj;
out.dprime  = dprime;
out.cvAcc   = cvAcc;
out.motifRank = motifRank;
out.time    = Htime(:)';  % 1 x T copy
out.params  = prm;

% ---------------- plotting ----------------
if prm.DoPlots
    fb = @(x, y1, y2, a, c) patch([x(:)' fliplr(x(:)')], ...
                                  [y1(:)' fliplr(y2(:)')], ...
                                  c, 'EdgeColor','none', 'FaceAlpha',a);

    figVis = 'on';
    if ~prm.visibleFigs, figVis = 'off'; end

    % Projections: Go vs NoGo
    h_proj = figure('Name','Coding projection (Go vs NoGo)', 'Visible', figVis); hold on;
    m_go   = mean(proj(goMask,:),  1, 'omitnan');
    se_go  = std( proj(goMask,:),  0, 1, 'omitnan') / sqrt(sum(goMask));
    m_ng   = mean(proj(nogoMask,:),1, 'omitnan');
    se_ng  = std( proj(nogoMask,:),0, 1, 'omitnan') / sqrt(sum(nogoMask));

    col_go = [0.8 0   0];   % red
    col_ng = [0   0   0.8]; % blue

    fb(Htime, m_go-se_go, m_go+se_go, 0.3, col_go);
    plot(Htime, m_go, 'Color', col_go, 'LineWidth', 2);

    fb(Htime, m_ng-se_ng, m_ng+se_ng, 0.3, col_ng);
    plot(Htime, m_ng, 'Color', col_ng, 'LineWidth', 2);

    xline(0,'k:'); xline(2,'k:');
    xlabel('Time (s)'); ylabel('Projection');
    ylim(prm.projCdYlim);
    legend({'Go ± SE','Go','NoGo ± SE','NoGo'});
    box on;

    % d' vs time
    h_dprime = figure('Name','Time-resolved d'' (per bin)', 'Visible', figVis);
    plot(Htime, dprime, 'LineWidth',2);
    xline(0,'k:'); xline(2,'k:');
    ylim(prm.dPrimeYlim);
    xlabel('Time (s)'); ylabel('d''');
    box on;

    % Wz heatmap
    h_Wz = figure('Name','Motif weights heatmap (standardized)', 'Visible', figVis);
    imagesc(Htime, 1:K, Wz); axis xy;
    xlabel('Time (s)'); ylabel('Motif #');
    colorbar;
    clim(prm.cdWeightCscale);
    title('Coding weights (standardized units; + = Go)');
    hold on;
    ylims = ylim;
    plot([0 0; 2 2],[ylims; ylims],'k:');

    % print
    if prm.printPlots
        if exist(figSaveDir,'dir')~=7, mkdir(figSaveDir); end
        header = extract_date_animalID_header(figSaveDir);
        timestampStr = datestr(now, 'mmddyy_HHMMSS');

        figSaveName_proj   = sprintf('proj_score_CD_bins_%s_%s', header, timestampStr);
        figSaveName_dprime = sprintf('dPrime_CD_bins_%s_%s',     header, timestampStr);
        figSaveName_Wz     = sprintf('wZ_CD_bins_%s_%s',         header, timestampStr);

        print(h_proj,   fullfile(figSaveDir, figSaveName_proj),  '-dpdf','-painters','-bestfit');
        print(h_dprime, fullfile(figSaveDir, figSaveName_dprime),'-dpdf','-painters','-bestfit');
        print(h_Wz,     fullfile(figSaveDir, figSaveName_Wz),    '-dpdf','-painters','-bestfit');
    end

    if ~prm.visibleFigs
        close(h_proj); close(h_dprime); close(h_Wz);
    end
end

end
