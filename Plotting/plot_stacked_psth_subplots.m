function hFig = plot_stacked_psth_subplots(H, t, varargin)
%PLOT_STACKED_PSTH_SUBPLOTS  Stack K motif time series in Kx1 subplots (top-down).
%
% hFig = PLOT_STACKED_PSTH_SUBPLOTS(H, t, 'CueLines', [0 2], 'LinkX', true, 'FigVisible','on')
%
% INPUTS
%   H           K x T matrix (motif x time bins)
%   t           1 x T time vector (optional). If empty or omitted, uses 1:T.
%
% NAME-VALUE OPTIONS
%   'CueLines'  : vector of x positions to draw as dotted lines (default: [])
%   'LinkX'     : logical, link x-axes across subplots (default: true)
%   'FigVisible': 'on' or 'off' (default: 'on')
%   'XLim'      : [xmin xmax] for all panels (default: [])
%   'YLabel'    : char/string for each panel’s ylabel prefix (default: 'Motif')
%
% OUTPUT
%   hFig        : figure handle
%
% Notes
%   - Subplot order is top-down: row 1 at the top, row K at the bottom.
%   - Each subplot auto-scales its own y-limits (no normalization).
%
% Junchol-ready utility, 2025.

% ---------- parse inputs ----------
if nargin < 2 || isempty(t)
    t = 1:size(H,2);
end
p = inputParser;
p.addParameter('CueLines', [], @(v)isnumeric(v)||isempty(v));
p.addParameter('LinkX', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('FigVisible', 'on', @(s)ischar(s)||isstring(s));
p.addParameter('XLim', [], @(v)isnumeric(v)&&numel(v)==2 || isempty(v));
p.addParameter('YLabel', 'Motif', @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
prm = p.Results;

% ---------- checks ----------
[K, T] = size(H);
assert(numel(t)==T, 'Length of t must match size(H,2).');

% ---------- figure ----------
hFig = figure('Name','Stacked PSTHs','Color','w','Visible',char(prm.FigVisible));

% ---------- plot, top-down ----------
ax = gobjects(K,1);
for k = 1:K
    ax(k) = subplot(K,1,k);
    plot(t, H(k,:), 'k', 'LineWidth', 1); hold on;

    % cue lines
    if ~isempty(prm.CueLines)
        for xc = prm.CueLines(:)'
            xline(xc, 'k:', 'LineWidth', 1);
        end
    end

    % labels and cosmetics
    %ylabel(sprintf('%s %d', string(prm.YLabel), k), 'Interpreter','none');
    set(ax(k), 'YTickLabel', []); % mute numeric y-ticks
    box off;

    % x limits
    if ~isempty(prm.XLim)
        xlim(ax(k), prm.XLim);
    end

    % Only bottom subplot gets x-label
    if k ~= K
        set(ax(k), 'XTickLabel', []);
    else
        xlabel('Time (s)', 'Interpreter','none');
    end
end

% link x if requested
if prm.LinkX && K > 1
    linkaxes(ax, 'x');
end

drawnow;
end
