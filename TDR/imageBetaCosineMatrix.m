function hFig = imageBetaCosineMatrix(S_motif, sessLabels, varargin)
% imageBetaCosineMatrix  Visualize β cosine-similarity across sessions.
%
% hFig = imageBetaCosineMatrix(S_motif, sessLabels, ...
%           'motifIdx', [], ...
%           'CLim', [0.5 1], ...
%           'Colormap', parula, ...
%           'Parent', [], ...
%           'FigSaveBaseDir', '')
%
% INPUTS
%   S_motif    : 1xK cell, each cell is an SxS cosine-similarity matrix (double, 0..1)
%   sessLabels : 1xS cellstr of session labels (used for both x/y ticks)
%
% NAME–VALUE OPTIONS
%   'motifIdx'       : [] (default) | scalar index 1..K. If empty, plot all motifs.
%   'CLim'           : [min max] for color axis (default [0.5 1])
%   'Colormap'       : colormap matrix | name | function handle (default parula)
%   'Parent'         : axes handle to plot into. If empty, new figure(s) are created.
%   'FigSaveBaseDir' : char/str; if nonempty, saves PDF(s) into <dir>/<animalID>/
%
% OUTPUT
%   hFig : handle to the figure (if plotting one motif or using Parent),
%          or a vector of figure handles if plotting all motifs with new figures.

% ---- Parse inputs ----
p = inputParser;
p.addParameter('motifIdx', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x)));
p.addParameter('CLim', [0.5 1], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('Colormap', parula, @(x) isnumeric(x) || isa(x,'function_handle') || ischar(x) || isstring(x));
p.addParameter('Parent', [], @(x) isempty(x) || ishghandle(x,'axes'));
p.addParameter('FigSaveBaseDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

% ---- Basic checks ----
assert(iscell(S_motif) && ~isempty(S_motif), 'S_motif must be a non-empty 1xK cell.');
K = numel(S_motif);
assert(iscellstr(sessLabels) && ~isempty(sessLabels), 'sessLabels must be a non-empty cellstr.');
S = numel(sessLabels);

% motif selection
if ~isempty(opt.motifIdx)
    assert(opt.motifIdx>=1 && opt.motifIdx<=K, 'motifIdx out of range (1..K).');
    motifList = opt.motifIdx(:).';
else
    motifList = 1:K;
end

% extract an animal/subject ID from labels, if possible
mIdTok = regexp(sessLabels{1}, '(m\d{3,5})', 'tokens', 'once');
if ~isempty(mIdTok)
    mId = mIdTok{1};
else
    mId = 'animal';
end

% Pre-allocate outputs
hFigs = gobjects(1, numel(motifList));

for ii = 1:numel(motifList)
    k = motifList(ii);

    Sk = S_motif{k};
    assert(ismatrix(Sk) && size(Sk,1)==S && size(Sk,2)==S, ...
        'S_motif{%d} must be SxS where S = numel(sessLabels).', k);

    % ---- Title string ----
    titleStr = sprintf('beta cosine similarity across sessions — motif %d of %s', k, mId);

    % ---- Axes / Figure ----
    if isempty(opt.Parent)
        hFigs(ii) = figure('Color','w','Visible','off'); %#ok<LAXES>
        ax = axes('Parent', hFigs(ii));
    else
        ax = opt.Parent;
        hFigs(ii) = ancestor(ax, 'figure');
    end

    % ---- Image ----
    imagesc(ax, Sk);
    axis(ax, 'image');
    set(ax, 'TickDir','out', 'Box','off');
    colormap(ax, opt.Colormap);
    caxis(ax, opt.CLim);
    colorbar(ax);

    % ---- Labels (turn off interpreters) ----
    set(ax, 'XTick', 1:S, 'YTick', 1:S, ...
        'XTickLabel', sessLabels, 'YTickLabel', sessLabels, ...
        'XTickLabelRotation', 45, ...
        'TickLabelInterpreter','none');

    title(ax, titleStr, 'Interpreter','none');

    % ---- Optional save ----
    if ~isempty(opt.FigSaveBaseDir)
        figSaveDir = fullfile(char(opt.FigSaveBaseDir), mId);
        if exist(figSaveDir, 'dir') ~= 7
            mkdir(figSaveDir);
        end
        baseName = sprintf('%s_beta_cosine_similarity_across_sessions_motif_%02d', mId, k);
        outPDF   = fullfile(figSaveDir, [baseName '.pdf']);
        print(hFigs(ii), outPDF, '-dpdf', '-painters', '-bestfit');
    end
end

% output: single handle or vector
if numel(hFigs)==1
    hFig = hFigs;
else
    hFig = hFigs;
end
end
