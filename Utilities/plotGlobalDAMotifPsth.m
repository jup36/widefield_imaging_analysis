function out = plotGlobalDAMotifPsth(S, motifK, varargin)
%PLOTGLOBALDAMOTIFPSTH
%   Overlay one motif's H PSTH with the global DA PSTH, on the time axis
%   the cross-correlogram was actually computed on.
%
%   The PSTHs are stored by globalDA_motifH_perTrial_xcorr_func whether or
%   not PSTH subtraction was applied, so this works either way -- and when
%   subtraction IS on, these are exactly the traces that were removed
%   before correlating. Reading the correlogram against them is the quickest
%   way to see whether a peak reflects residual co-fluctuation or leftover
%   shared task locking.
%
%   BOTH SIGNALS ARE Z-SCORED (globally, inside stack_trials_H), so they
%   share units and can go on one axis. If that ever changes, or if you
%   want raw dF/F, set 'twoAxes' true for a right-hand axis instead of
%   silently plotting different units on the same scale.
%
%   out = plotGlobalDAMotifPsth(S, motifK, ...)
%
% INPUTS
%   S      : output of globalDA_motifH_perTrial_xcorr_func
%   motifK : motif row to plot
%
% NAME-VALUE
%   'stream'    : 'hit' (default) | 'cr'
%   'twoAxes'   : false (default). true -> DA on a right-hand y-axis.
%   'hColor'    : [0.00 0.45 0.74]  (blue)
%   'daColor'   : [0.47 0.25 0.80]  (purple -- distinct from the trial-type
%                 palette, since this is not a trial type)
%   'lineWidth' : 2
%   'ax'        : axes to draw into; default new figure
%   'titleStr'  : '' (auto)
%
% OUTPUT
%   .timeX, .psthH, .psthDA, .fig, .ax

p = inputParser;
p.addParameter('stream', 'hit', @(s) any(strcmpi(string(s), ["hit","cr"])));
p.addParameter('twoAxes', false, @(x) islogical(x) && isscalar(x));
p.addParameter('hColor', [0.00 0.45 0.74], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('daColor', [0.47 0.25 0.80], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('lineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
p.addParameter('titleStr', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

stream = lower(char(string(opt.stream)));

%% ---- pull ----
fH  = ['psthH_'  stream];
fDA = ['psthDA_' stream];
assert(isfield(S.obs, fH) && isfield(S.obs, fDA), ...
    'S.obs.%s / %s not found -- was the "%s" stream computed?', fH, fDA, stream);

psthHall = S.obs.(fH);
K = size(psthHall, 1);
assert(motifK >= 1 && motifK <= K, 'motifK must be in 1..%d.', K);

psthH  = psthHall(motifK, :);
psthDA = S.obs.(fDA);

% The x-axis is the stack_trials_H window centres -- the same grid the
% correlogram lags are defined on. Falling back to bin index would make the
% two figures silently incomparable, so this errors instead.
assert(isfield(S.meta, 'winCtrs') && numel(S.meta.winCtrs) == numel(psthH), ...
    'S.meta.winCtrs missing or wrong length -- cannot build a time axis.');
timeX = double(S.meta.winCtrs(:)');

nTr = NaN;
if isfield(S.meta, ['nTrials_' stream]), nTr = S.meta.(['nTrials_' stream]); end

%% ---- plot ----
if isempty(opt.ax)
    fig = figure('Color', 'w');
    ax  = axes(fig);
else
    ax  = opt.ax;
    fig = ancestor(ax, 'figure');
end

if opt.twoAxes
    yyaxis(ax, 'left');
    hH = plot(ax, timeX, psthH, '-', 'Color', opt.hColor, 'LineWidth', opt.lineWidth);
    ylabel(ax, sprintf('Motif %d H (z)', motifK));
    ax.YColor = opt.hColor;

    yyaxis(ax, 'right');
    hDA = plot(ax, timeX, psthDA, '-', 'Color', opt.daColor, 'LineWidth', opt.lineWidth);
    ylabel(ax, 'Global DA (z)');
    ax.YColor = opt.daColor;

    yyaxis(ax, 'left');
else
    hold(ax, 'on');
    hH  = plot(ax, timeX, psthH,  '-', 'Color', opt.hColor,  'LineWidth', opt.lineWidth);
    hDA = plot(ax, timeX, psthDA, '-', 'Color', opt.daColor, 'LineWidth', opt.lineWidth);
    ylabel(ax, 'PSTH (z-scored)');
    yline(ax, 0, '-', 'Color', [0.82 0.82 0.82], 'HandleVisibility', 'off');
end

xline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');

xlabel(ax, 'Time from tone onset (s)');
legend(ax, [hH hDA], {sprintf('Motif %d H', motifK), 'Global DA'}, 'Location', 'best');

if strlength(string(opt.titleStr)) > 0
    title(ax, char(opt.titleStr), 'Interpreter', 'none');
else
    psthTag = 'PSTH-subtracted run';
    if isfield(S.params, 'doPSTHSubtraction') && ~S.params.doPSTHSubtraction
        psthTag = 'no PSTH subtraction';
    end
    title(ax, sprintf('Motif %d H vs. global DA | %s trials (n=%d) | %s', ...
        motifK, upper(stream), nTr, psthTag), 'Interpreter', 'none');
end

xlim(ax, [timeX(1) timeX(end)]);
set(ax, 'TickDir', 'out');
grid(ax, 'on');
box(ax, 'off');
hold(ax, 'off');

out = struct('timeX', timeX, 'psthH', psthH, 'psthDA', psthDA, ...
    'motifK', motifK, 'stream', stream, 'fig', fig, 'ax', ax);
end