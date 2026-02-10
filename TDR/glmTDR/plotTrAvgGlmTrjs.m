function plotTrAvgGlmTrjs(zTrj, trIC, targetColumns, timestamps, varargin)
%PLOTTRAVGGLMTRJS  Plot trial-averaged GLM-TDR projected trajectories.
%
% MODS (2026-02-09)
%   - Marks trajectory start with OPEN circle, end with FILLED circle.
%   - Adds a time-flow cue by drawing the trajectory as short segments whose
%     color gradually ramps from faint -> solid.
%
% NOTE ON "ALPHA"
%   MATLAB line objects don't support per-vertex alpha. The most robust way
%   is to draw the trajectory as many short line segments and set each
%   segment's Color to an RGB that linearly interpolates from white->baseColor.
%   (So it LOOKS faint at start and solid at end.)
%
%   If you truly want transparency alpha, you can instead use patch objects,
%   but that’s more brittle across versions. This implementation is stable.

% -------------------- parse --------------------
p = inputParser;
p.addParameter('tBounds', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1)<x(2)));
p.addParameter('colorMapC', {}, @(c) iscell(c) || isstring(c));
p.addParameter('smoothingFactor', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('LineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.addParameter('TitleStr', '', @(s) ischar(s) || isstring(s));
p.addParameter('LegendC', {}, @(c) isempty(c) || iscell(c) || isstring(c));

% new
p.addParameter('StartMarkerSize', 48, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('EndMarkerSize',   48, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('FadeToWhite',      0.85, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1); % 0=no fade, 1=very faint start
p.addParameter('FadeN',            80, @(x) isnumeric(x) && isscalar(x) && x>=5);           % number of segments used for fade
p.parse(varargin{:});
opt = p.Results;

% -------------------- sanity --------------------
assert(ndims(zTrj)==3, 'zTrj must be N x T x D.');
[N, T, D] = size(zTrj);

timestamps = timestamps(:)'; % 1 x T
assert(numel(timestamps)==T, 'timestamps must have length T matching size(zTrj,2).');

if ~iscell(trIC), error('trIC must be a cell array of trial-index logical vectors.'); end
C = numel(trIC);

targetColumns = targetColumns(:)';
assert(all(targetColumns>=1 & targetColumns<=D), 'targetColumns out of range for zTrj third dimension.');

nDimPlot = numel(targetColumns);
if nDimPlot >= 3
    useCols = targetColumns(1:3);
elseif nDimPlot == 2
    useCols = targetColumns(1:2);
elseif nDimPlot == 1
    useCols = targetColumns(1);
else
    error('targetColumns is empty.');
end

% -------------------- time bounds --------------------
tMask = true(1,T);
if ~isempty(opt.tBounds)
    tMask = (timestamps >= opt.tBounds(1)) & (timestamps <= opt.tBounds(2));
end
tIdx = find(tMask);
if isempty(tIdx)
    warning('tBounds excluded all timestamps. Nothing to plot.');
    return;
end

% -------------------- colors --------------------
if isempty(opt.colorMapC)
    cmap = lines(max(C,1));
    colorMapC = arrayfun(@(i) cmap(i,:), 1:C, 'UniformOutput', false);
else
    colorMapC = opt.colorMapC;
    if isstring(colorMapC), colorMapC = cellstr(colorMapC); end
    if numel(colorMapC) < C
        last = colorMapC{end};
        for i = (numel(colorMapC)+1):C
            colorMapC{i} = last; %#ok<AGROW>
        end
    end
end

% -------------------- legend --------------------
if isempty(opt.LegendC)
    legC = arrayfun(@(i) sprintf('cond%d', i), 1:C, 'UniformOutput', false);
else
    legC = cellstr(string(opt.LegendC));
    if numel(legC) < C
        for i = (numel(legC)+1):C
            legC{i} = sprintf('cond%d', i); %#ok<AGROW>
        end
    end
end

% -------------------- figure --------------------
if opt.MakeFigure
    figure('Color','w'); hold on;
end

% handles for legend (one per condition)
hLeg = gobjects(1,C);

% -------------------- main loop: trial-avg trajectories --------------------
for c = 1:C
    trI = trIC{c};
    if isempty(trI), continue; end
    if ~islogical(trI)
        tmp = false(N,1);
        tmp(trI(:)) = true;
        trI = tmp;
    end
    assert(numel(trI)==N, 'Each trIC{%d} must be logical [N x 1].', c);
    if ~any(trI), continue; end

    % Extract [nTrials x T x nDim]
    Zc = zTrj(trI, :, useCols);

    % Trial-average => [T x nDim] (omit NaNs)
    Zm = squeeze(mean(Zc, 1, 'omitnan'));  % [T x nDim] or [T x 1]
    if isvector(Zm), Zm = Zm(:); end

    % Restrict time
    Zm = Zm(tIdx, :);
    tPlot = timestamps(tIdx);

    % Smooth along time
    if opt.smoothingFactor > 0
        Zm = smooth2a(Zm, opt.smoothingFactor, 0);
    end

    baseCol = colorMapC{c};

    % ---- draw faded trajectory as segments (faint -> solid) ----
    % choose points to use for segmenting (downsample in time if long)
    nPts = size(Zm,1);
    if nPts < 2, continue; end

    % Use at most opt.FadeN segments for speed/clarity
    nSeg = min(opt.FadeN, nPts-1);

    % indices spanning 1..nPts
    idxPts = unique(round(linspace(1, nPts, nSeg+1)));
    if numel(idxPts) < 2, idxPts = 1:nPts; end

    % fade parameter: start near white, end at base color
    % w=0 -> solid base; w=1 -> white
    w0 = opt.FadeToWhite;
    wVec = linspace(w0, 0, numel(idxPts)-1);

    % plot segments
    for s = 1:(numel(idxPts)-1)
        i1 = idxPts(s);
        i2 = idxPts(s+1);

        segCol = blend_to_white(baseCol, wVec(s));

        if nDimPlot == 1
            hh = plot(tPlot([i1 i2]), Zm([i1 i2],1), ...
                'LineWidth', opt.LineWidth, 'Color', segCol);
        elseif nDimPlot == 2
            hh = plot(Zm([i1 i2],1), Zm([i1 i2],2), ...
                'LineWidth', opt.LineWidth, 'Color', segCol);
        else
            hh = plot3(Zm([i1 i2],1), Zm([i1 i2],2), Zm([i1 i2],3), ...
                'LineWidth', opt.LineWidth, 'Color', segCol);
        end

        % keep one handle per condition for legend (use the LAST segment so it's solid)
        if s == (numel(idxPts)-1)
            hLeg(c) = hh;
        end
    end

    % ---- start/end markers ----
    pStart = Zm(1,:);
    pEnd   = Zm(end,:);

    if nDimPlot == 1
        % start: open circle
        scatter(tPlot(1),  pStart(1), opt.StartMarkerSize, ...
            'MarkerEdgeColor', baseCol, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        % end: filled circle
        scatter(tPlot(end), pEnd(1),   opt.EndMarkerSize, ...
            'MarkerEdgeColor', baseCol, 'MarkerFaceColor', baseCol, 'LineWidth', 1.0);

        xlabel('time');
        ylabel(sprintf('proj (axis %d)', useCols(1)));

    elseif nDimPlot == 2
        scatter(pStart(1), pStart(2), opt.StartMarkerSize, ...
            'MarkerEdgeColor', baseCol, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        scatter(pEnd(1),   pEnd(2),   opt.EndMarkerSize, ...
            'MarkerEdgeColor', baseCol, 'MarkerFaceColor', baseCol, 'LineWidth', 1.0);

        xlabel(sprintf('axis %d', useCols(1)));
        ylabel(sprintf('axis %d', useCols(2)));
        axis tight;

    else
        scatter3(pStart(1), pStart(2), pStart(3), opt.StartMarkerSize, ...
            'MarkerEdgeColor', baseCol, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        scatter3(pEnd(1),   pEnd(2),   pEnd(3),   opt.EndMarkerSize, ...
            'MarkerEdgeColor', baseCol, 'MarkerFaceColor', baseCol, 'LineWidth', 1.0);

        xlabel(sprintf('axis %d', useCols(1)));
        ylabel(sprintf('axis %d', useCols(2)));
        zlabel(sprintf('axis %d', useCols(3)));
        grid on;
        axis vis3d;
    end
end

% -------------------- cosmetics --------------------
if ~isempty(opt.TitleStr)
    title(opt.TitleStr, 'Interpreter','none');
end

hasLeg = isgraphics(hLeg);
if any(hasLeg)
    legend(hLeg(hasLeg), legC(hasLeg), 'Location','best', 'Interpreter','none');
end

end

%% ---------------- HELPER ----------------
function colOut = blend_to_white(colIn, w)
%BLEND_TO_WHITE  Return (1-w)*colIn + w*[1 1 1], clamped to [0,1]
colIn = colIn(:)';
if numel(colIn) ~= 3, colIn = [0 0 0]; end
w = max(min(w,1),0);
colOut = (1-w)*colIn + w*[1 1 1];
colOut = max(min(colOut,1),0);
end
