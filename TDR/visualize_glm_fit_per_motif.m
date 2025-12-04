function Yhat = visualize_glm_fit_per_motif(model, Xz, YbigVal, varargin)
% VISUALIZE_GLM_FIT_PER_MOTIF
% Reconstructs motif-wise predictions and visualizes whole-trace and
% best-matching snippet overlays.
%
% Usage
%   Yhat = visualize_glm_fit_per_motif(model, Xz, YbigVal, ...
%              'Time', t, 'WinSec', 5, 'Fs', [], 'Motifs', [], ...
%              'Link', 'identity', 'SaveFigDir', '')
%
% Inputs
%   model   : either struct with fields
%               .beta (P x K), optional .bias (1 x K or K x 1),
%               optional .invlink (function handle)
%             OR a numeric beta (P x K). If numeric, you can pass 'Bias', b0.
%   Xz      : T x P design matrix (z-scored predictors)
%   YbigVal : T x K response matrix (one column per motif)
%
% Name–Value options
%   'Bias'      : 1xK bias (intercept); default zeros if absent
%   'InvLink'   : function handle for inverse link (default @identity)
%   'Link'      : 'identity'|'exp'|'logistic' (overridden by InvLink if given)
%   'Time'      : T x 1 vector of timestamps (sec). If empty, uses samples.
%   'Fs'        : sampling rate (Hz). Used only if 'Time' empty & WinSec given.
%   'WinSec'    : snippet window length in seconds (or samples if no time)
%                 default: 5 sec (or 1000 samples if no time/Fs)
%   'Motifs'    : vector of motif indices to visualize (default: all)
%   'Smooth'    : optional smoothing (movmean window in samples). default: 0 (off)
%   'SaveFigDir': folder to save figures ('' = don’t save)
%
% Output
%   Yhat : T x K matrix of reconstructed predictions

% ---------- parse inputs
p = inputParser;
p.addParameter('Bias', [], @(x)isnumeric(x));
p.addParameter('InvLink', [], @(f)isempty(f)||isa(f,'function_handle'));
p.addParameter('Link', 'identity', @(s)any(strcmpi(s,{'identity','exp','logistic'})));
p.addParameter('Time', [], @(x)isempty(x)||isvector(x));
p.addParameter('Fs', [], @(x)isempty(x)||isscalar(x));
p.addParameter('WinSec', [], @(x)isempty(x)||isscalar(x));
p.addParameter('Motifs', [], @(x)isempty(x)||isvector(x));
p.addParameter('Smooth', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('SaveFigDir','',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

% ---------- unpack model
if isstruct(model)
    beta = model.beta;
    if isfield(model,'bias'), b0 = model.bias(:)'; else, b0 = []; end
    if isempty(opt.InvLink) && isfield(model,'invlink'), invlink = model.invlink; else, invlink = []; end
else
    beta = model;
    b0   = opt.Bias;
    invlink = opt.InvLink;
end
[P,K] = size(beta);
if isempty(b0), b0 = zeros(1,K); end
b0 = b0(:)';           % row

% pick inverse link
if isempty(invlink)
    switch lower(opt.Link)
        case 'identity', invlink = @(x)x;
        case 'exp',      invlink = @(x)exp(x);
        case 'logistic', invlink = @(x)1./(1+exp(-x));
    end
end

[T,Px] = size(Xz);  %#ok<ASGLU>
assert(Px==P,'Size mismatch: Xz is T×%d but beta is %d×K',Px,P);

assert(size(YbigVal,1)==T,'YbigVal must have T rows.');

% which motifs
if isempty(opt.Motifs), kmask = 1:K; else, kmask = opt.Motifs(:)'; end

% time base / window
t = opt.Time;
if isempty(t)
    if ~isempty(opt.Fs)
        t = (0:T-1)'/opt.Fs;
    else
        t = (1:T)'; % samples
    end
end

if isempty(opt.WinSec)
    % default window = 5s if we have seconds, else ~1000 samples
    if max(t) > T   % assume t is in seconds (not samples)
        win = 5;   % seconds
    else
        win = 1000; % samples
    end
else
    win = opt.WinSec;
end

% convert win to samples if t is seconds
if max(t) > T
    % seconds scale → samples via median dt
    dt  = median(diff(t),'omitnan');
    winSamp = max(5, round(win/dt));
else
    winSamp = max(5, round(win));
end

% optional smoothing window
smw = round(opt.Smooth);

% ---------- predict
linpred = Xz*beta + repmat(b0, T, 1);  % T×K
Yhat    = invlink(linpred);            % T×K

% optional smoothing (light) for viewing only
if smw>1
    Yplot  = movmean(YbigVal, smw, 1, 'omitnan');
    Yhplot = movmean(Yhat,    smw, 1, 'omitnan');
else
    Yplot  = YbigVal;
    Yhplot = Yhat;
end

% ---------- visualization per motif
for k = kmask
    y   = Yplot(:,k);
    yh  = Yhplot(:,k);

    % metrics
    c   = corr(y, yh, 'rows','pairwise');
    ssr = nansum((y - yh).^2);
    sst = nansum((y - nanmean(y)).^2);
    R2  = 1 - ssr/sst;

    % moving correlation to choose best snippet
    mc  = movcorr(y, yh, winSamp, 'Endpoints','shrink');  % T×1
    [~, iBest] = max(mc);
    i1 = max(1, iBest - floor(winSamp/2));
    i2 = min(T, i1 + winSamp - 1);
    i1 = max(1, i2 - winSamp + 1); % enforce exact length when possible

    % ---- figure
    fig = figure('Color','w','Name',sprintf('Motif %d fit',k));
    tiledlayout(fig, 2, 1, 'Padding','compact', 'TileSpacing','compact');

    % (1) Full trace
    ax1 = nexttile;
    plot(ax1, t, y,  'Color',[0.2 0.4 0.8 0.8], 'LineWidth', 1); hold(ax1,'on');
    plot(ax1, t, yh, 'Color',[0.8 0.2 0.2 0.8], 'LineWidth', 1);
    xlabel(ax1, 'Time'); ylabel(ax1, sprintf('Motif %d',k));
    title(ax1, sprintf('Whole trace — corr = %.3f, R^2 = %.3f', c, R2));
    legend(ax1, {'Y','Ŷ'}, 'Location','best'); box(ax1,'off'); grid(ax1,'on');

    % (2) Best snippet
    ax2 = nexttile;
    plot(ax2, t(i1:i2), y(i1:i2),  'Color',[0.2 0.4 0.8 0.9], 'LineWidth', 1.5); hold(ax2,'on');
    plot(ax2, t(i1:i2), yh(i1:i2), 'Color',[0.8 0.2 0.2 0.9], 'LineWidth', 1.5);
    xlabel(ax2, 'Time'); ylabel(ax2, sprintf('Motif %d',k));
    title(ax2, sprintf('Best %d-sample window (max mov. corr = %.3f)', i2-i1+1, max(mc,'omitnan')));
    legend(ax2, {'Y','Ŷ'}, 'Location','best'); box(ax2,'off'); grid(ax2,'on');

    % save if requested
    if ~isempty(opt.SaveFigDir)
        if ~exist(opt.SaveFigDir,'dir'), mkdir(opt.SaveFigDir); end
        fn = fullfile(opt.SaveFigDir, sprintf('glm_fit_motif_%02d.png',k));
        exportgraphics(fig, fn, 'Resolution', 200);
    end
end
end


function r = movcorr(x, y, w, varargin)
%MOVCORR  Moving window correlation (compatible with pre-R2023b MATLAB)
%
%   r = movcorr(x, y, w)
%
%   Computes Pearson correlation between x and y in a sliding window
%   of width w (in samples).  Returns vector of same length as x.
%
%   Optional 'Endpoints' behavior:
%       'shrink'  – compute correlation using smaller windows near edges (default)
%       'discard' – return NaN near edges
%
%   (Lightweight stand-in for MATLAB's built-in movcorr.)

    if nargin < 3
        error('movcorr(x, y, w) requires three inputs');
    end
    if isempty(x) || isempty(y)
        r = NaN(size(x));
        return;
    end

    x = x(:); y = y(:);
    N = numel(x);
    r = NaN(N,1);

    % parse optional 'Endpoints'
    endpoints = 'shrink';
    if nargin > 3 && ischar(varargin{1})
        endpoints = lower(varargin{1});
    end

    halfw = floor(w/2);
    for i = 1:N
        i1 = max(1, i - halfw);
        i2 = min(N, i + halfw);
        xi = x(i1:i2);
        yi = y(i1:i2);
        if numel(xi) > 2 && all(isfinite(xi)) && all(isfinite(yi))
            r(i) = corr(xi, yi, 'rows', 'pairwise');
        elseif strcmp(endpoints,'discard')
            r(i) = NaN;
        end
    end
end
