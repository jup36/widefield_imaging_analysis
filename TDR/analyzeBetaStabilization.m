function out = analyzeBetaStabilization(S_motif, varargin)
% ANALYZEBETASTABILIZATION  Quantify when β–cosine stabilizes (one or all motifs).
%
% out = analyzeBetaStabilization(S_motif, ...
%       'motifIdx', [], ...                 % [] => all motifs, else scalar 1..K
%       'sessionLabels', {...}, ...         % 1xT cellstr (default {'s1',...})
%       'lateN', 3, ...                     % #final sessions defining plateau
%       'threshMode','rel', ...             % 'rel' (default) | 'abs'
%       'thresh', 0.8, ...                  % rel: fraction in [0,1]; abs: value
%       'smoothSigma', 1, ...               % display smoothing sigma (0 => off)
%       'onsetSmoothSigma', 0, ...          % smoothing used for onset detection
%       'onsetConsec', 1, ...               % require N consecutive ≥ threshold
%       'onsetStart', 1, ...                % first session index to consider
%       'plot', true, ...                   % make diagnostic plots (invisible)
%       'figSaveBaseDir', '', ...           % if nonempty, save PDFs to this dir
%       'showColorbar', false)              % show colorbar only in diagnostic fig
%
% INPUT
%   S_motif : 1xK cell, each K cell is TxT cosine similarity; or a single TxT.
%
% OUTPUT (struct per analyzed motif)
%   .motifIdx, .v, .v_smooth
%   .onset_idx, .onset_label
%   .plateau_mean, .base_mean, .thresh_value
%   .fig (diagnostic figure, invisible), .figFile (saved PDF path)

% -------- Parse inputs --------
p = inputParser;
p.addParameter('motifIdx', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x)));
p.addParameter('sessionLabels', {}, @(x) iscell(x) || isempty(x));
p.addParameter('lateN', 3, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('threshMode', 'rel', @(s) any(strcmpi(s,{'rel','abs'})));
p.addParameter('thresh', 0.8, @(x) isnumeric(x) && isscalar(x));
p.addParameter('smoothSigma', 1, @(x) isnumeric(x) && isscalar(x) && x>=0);
% Onset-detector knobs
p.addParameter('onsetSmoothSigma', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('onsetConsec', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('onsetStart', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
% plotting / saving
p.addParameter('plot', true, @(x) islogical(x) || ismember(x,[0 1]));
p.addParameter('figSaveBaseDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('showColorbar', false, @(x) islogical(x) || ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

% -------- Normalize S_motif to a cell array --------
if ~iscell(S_motif), S_motif = {S_motif}; end
assert(~isempty(S_motif), 'S_motif must be non-empty.');
K = numel(S_motif);

% -------- Basic size checks & labels --------
S1 = S_motif{1};
assert(ismatrix(S1) && size(S1,1)==size(S1,2), 'Each S must be square TxT.');
T = size(S1,1);
if isempty(opt.sessionLabels)
    opt.sessionLabels = arrayfun(@(i) sprintf('s%d', i), 1:T, 'uni', 0);
end
assert(numel(opt.sessionLabels)==T, 'sessionLabels must have T items.');

% -------- Motif selection --------
if isempty(opt.motifIdx)
    motifList = 1:K;
else
    assert(opt.motifIdx>=1 && opt.motifIdx<=K, 'motifIdx out of range (1..K).');
    motifList = opt.motifIdx;
end

% -------- Extract an ID (e.g., m####) from session labels (if present) --------
mIdTok = regexp(opt.sessionLabels{1}, '(m\d{3,5})', 'tokens', 'once');
if ~isempty(mIdTok), mId = mIdTok{1}; else, mId = 'animal'; end

% -------- Analyze each motif --------
out = repmat(struct( ...
    'motifIdx', [], ...
    'v', [], 'v_smooth', [], ...
    'onset_idx', NaN, 'onset_label', '', ...
    'plateau_mean', NaN, 'base_mean', NaN, 'thresh_value', NaN, ...
    'fig', [], 'figFile', ''), 1, numel(motifList));

for ii = 1:numel(motifList)
    k = motifList(ii);
    S = double(S_motif{k});
    assert(all(size(S)==[T T]), 'S_motif{%d} must be %dx%d to match sessionLabels.', k, T, T);

    % ---- v(i): mean similarity to the final lateN sessions (RAW) ----
    lateN = min(opt.lateN, T);
    v = nan(1,T);
    for i = 1:T
        v(i) = mean(S(i, T-lateN+1:T), 'omitnan');
    end

    % ---- Edge-safe smoothing for display; pin endpoints to raw ----
    v_smooth = v;
    if opt.smoothSigma > 0
        win = max(3, 2*ceil(3*opt.smoothSigma)+1);         % ~6*sigma+1 (odd)
        v_smooth = gaussSmoothReplicate(v, win);
        v_smooth([1 end]) = v([1 end]);
    end

    % ---- Baseline & plateau (SMOOTHED for reporting) ----
    plateau_idx = (T-lateN+1):T;
    plateau_mean = mean(v_smooth(plateau_idx), 'omitnan');
    early_ref = min(2, T-1);
    base_mean = mean(v_smooth(1:early_ref), 'omitnan');

    % ---- Threshold from RAW anchors (robust) ----
    plateau_raw = mean(v(end-lateN+1:end), 'omitnan');
    base_raw    = mean(v(1:early_ref), 'omitnan');
    switch lower(opt.threshMode)
        case 'rel'
            thresh_value = base_raw + opt.thresh*(plateau_raw - base_raw);
        case 'abs'
            thresh_value = opt.thresh;
    end

    % ---- Onset detection on lightly-smoothed (or raw) series ----
    v_onset = v;
    if opt.onsetSmoothSigma > 0
        win_on = max(3, 2*ceil(3*opt.onsetSmoothSigma)+1);
        v_onset = gaussSmoothReplicate(v_onset, win_on);
        v_onset([1 end]) = v([1 end]);
    end
    onset_idx = NaN;
    above = v_onset >= thresh_value;
    startI = max(1, round(opt.onsetStart));
    for i = startI : T - opt.onsetConsec + 1
        if all(above(i : i + opt.onsetConsec - 1))
            onset_idx = i; break
        end
    end
    onset_label = ''; if ~isnan(onset_idx), onset_label = opt.sessionLabels{onset_idx}; end

    % ---- Pack output ----
    out(ii).motifIdx = k;
    out(ii).v = v;
    out(ii).v_smooth = v_smooth;
    out(ii).onset_idx = onset_idx;
    out(ii).onset_label = onset_label;
    out(ii).plateau_mean = plateau_mean;
    out(ii).base_mean = base_mean;
    out(ii).thresh_value = thresh_value;
    out(ii).fig = [];
    out(ii).figFile = '';

    % ---- Plot (invisible) and optionally save ----
    if opt.plot
        hFig = figure('Color','w', 'Visible','off');
        tl = tiledlayout(hFig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

        % Heatmap (diagnostic only; never exported)
        ax1 = nexttile(tl);
        imagesc(ax1, S); axis(ax1, 'image'); colormap(ax1, parula);
        if opt.showColorbar, try, colorbar(ax1); catch, end, end
        title(ax1, 'Cosine similarity (S)', 'Interpreter','none');
        xticks(ax1, 1:T); yticks(ax1, 1:T);
        xticklabels(ax1, opt.sessionLabels); yticklabels(ax1, opt.sessionLabels);
        set(ax1, 'TickLabelInterpreter','none'); xtickangle(ax1, 45);

        % Stability curve + threshold & onset
        ax2 = nexttile(tl); hold(ax2, 'on');
        plot(ax2, 1:T, v, 'o-', 'LineWidth', 1, 'DisplayName','v (raw)');
        plot(ax2, 1:T, v_smooth, '-', 'LineWidth', 2, 'DisplayName','v (smooth)');
        yline(ax2, thresh_value, ':', 'DisplayName','threshold');
        if ~isnan(onset_idx), xline(ax2, onset_idx, ':', 'DisplayName','onset'); end
        grid(ax2, 'on'); box(ax2, 'off');
        set(ax2, 'TickLabelInterpreter','none');
        xlabel(ax2, 'session', 'Interpreter','none');
        ylabel(ax2, 'mean similarity to later sessions', 'Interpreter','none');
        ttl = sprintf('Motif %d — Onset: %s', k, onset_label);
        title(ax2, ttl, 'Interpreter','none');
        legend(ax2, 'Location','best', 'Box','off', 'Interpreter','none');

        out(ii).fig = hFig;

        % ---- Minimal export figure (extra margin; no tiledlayout/colorbar) ----
        if ~isempty(opt.figSaveBaseDir)
            figDir  = fullfile(char(opt.figSaveBaseDir), mId);
            if exist(figDir, 'dir') ~= 7, mkdir(figDir); end
            baseName = sprintf('%s_beta_stabilization_motif_%02d', mId, k);
            outPDF   = fullfile(figDir, [baseName '.pdf']);

            hSave = figure('Visible','off', 'Color','w', 'Renderer','opengl', ...
                           'Units','pixels', 'Position',[100 100 1400 520]); % more margin

            % Left panel (margins enlarged)
            axL = axes('Parent', hSave, 'Units','normalized', ...
                       'Position',[0.08 0.15 0.38 0.76]);
            imagesc(axL, S); axis(axL, 'image'); colormap(axL, parula);
            title(axL, 'Cosine similarity (S)', 'Interpreter','none');
            xticks(axL, 1:T); yticks(axL, 1:T);
            xticklabels(axL, opt.sessionLabels); yticklabels(axL, opt.sessionLabels);
            set(axL, 'TickLabelInterpreter','none'); xtickangle(axL, 45);
            box(axL,'off'); set(axL,'TickDir','out');

            % Right panel (margins enlarged)
            axR = axes('Parent', hSave, 'Units','normalized', ...
                       'Position',[0.53 0.15 0.40 0.76]); hold(axR,'on');
            plot(axR, 1:T, v, 'o-', 'LineWidth', 1, 'DisplayName','v (raw)');
            plot(axR, 1:T, v_smooth, '-', 'LineWidth', 2, 'DisplayName','v (smooth)');
            yline(axR, thresh_value, ':', 'DisplayName','threshold');
            if ~isnan(onset_idx), xline(axR, onset_idx, ':', 'DisplayName','onset'); end
            grid(axR, 'on'); box(axR,'off');
            set(axR, 'TickLabelInterpreter','none', 'TickDir','out');
            xlabel(axR, 'session', 'Interpreter','none');
            ylabel(axR, 'mean similarity to later sessions', 'Interpreter','none');
            title(axR, sprintf('Motif %d — Onset: %s', k, onset_label), 'Interpreter','none');
            legend(axR, 'Location','best', 'Box','off', 'Interpreter','none');

            % Defensively ensure no colorbar and flush updates
            delete(findall(hSave, 'Type', 'ColorBar'));
            drawnow limitrate nocallbacks

            % Export (vector if possible), avoiding -painters
            try
                exportgraphics(hSave, outPDF, 'ContentType','vector', 'BackgroundColor','white');
            catch
                try
                    print(hSave, outPDF, '-dpdf', '-opengl', '-r300');
                catch
                    print(hSave, strrep(outPDF,'.pdf','.png'), '-dpng', '-opengl', '-r300');
                end
            end
            close(hSave);
            out(ii).figFile = outPDF;
        end
    end
end
end

% ---------- helpers ----------
function y = gaussSmoothReplicate(x, win)
% Replicate-padded Gaussian smoothing (toolbox-free), win must be odd.
win = max(3, win + mod(win+1,2));            % force odd >=3
sigma = (win-1)/6;                            % ~6*sigma+1 rule
g = exp(-((-(win-1)/2:(win-1)/2).^2)/(2*sigma^2)); g = g/sum(g);
k = (win-1)/2;
xpad = [repmat(x(1),1,k), x, repmat(x(end),1,k)];
y = conv(xpad, g, 'same');
y = y(1+k:end-k);
end
