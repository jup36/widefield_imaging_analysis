function out = perMouseAcrossSessionPrjScoreTrajectories(prj_glmA, headerC, timestamps, varargin)
%PERMOUSEACROSSSESSIONPRJSCORETRAJECTORIES  Plot 1D projected scores across sessions for one mouse.
%
% out = perMouseAcrossSessionPrjScoreTrajectories(prj_glmA, headerC, timestamps, 'Name', value, ...)
%
% REQUIRED INPUTS
%   prj_glmA    : output of projectGlmRezC_toAnchors(...)
%   headerC     : [J x S] cell of session headers (e.g., 'm1045_122424' or 'm1613_050725-1')
%   timestamps  : [1 x T] or [T x 1] time vector aligned to 2nd dim of each Z in prj_glmA.*.ZC
%x
% NAME-VALUE (core)
%   'mouseId'         : "m1045" (required)
%   'projType'        : "global" (default) | "perMouse"
%   'targetName'      : "NoGoToneOn_1" (required; must exist in names list)
%   'trIdC'           : cell [J x S], each entry is struct of trial selectors (optional)
%   'trialField'      : "crI" (default). If trIdC provided, uses trIdC{j,s}.(trialField) to pick trials.
%
% NAME-VALUE (time + filtering)
%   'tBounds'         : [] (default) or [tMin tMax]
%   'smoothingFactor' : 0 (default) or positive integer (e.g. 5)
%   'dateLaterThan'   : [] (default) or "MMDDYY" (include sessions >= this date)
%   'dateEarlierThan' : [] (default) or "MMDDYY" (include sessions <= this date)
%   'day4MarkC'       : [] (default) or cell/string array {mouseId, "MMDDYY"; ...}
%                      If provided, overrides dateLaterThan by using the mouse-specific Day4 date.
%                      If mouseId is not found in day4MarkC, throws an error.
%
% NAME-VALUE (visual)
%   'FadeToWhite'     : 0.85 (default). Early sessions are blended toward white by this amount.
%   'LineWidth'       : 2 (default)
%   'MakeFigure'      : true (default)
%
% NAME-VALUE (save)
%   'figSaveDir'      : '' (default). If non-empty, saves PDF to this directory.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('prj_glmA', @(s)isstruct(s));
p.addRequired('headerC', @(c) iscell(c));
p.addRequired('timestamps', @(t) isnumeric(t) && isvector(t));

p.addParameter('mouseId', "", @(s)ischar(s)||isstring(s));
p.addParameter('projType', "global", @(s)ischar(s)||isstring(s));
p.addParameter('targetName', "", @(s)ischar(s)||isstring(s));

p.addParameter('trIdC', {}, @(c) isempty(c) || iscell(c));
p.addParameter('trialField', "", @(s)ischar(s)||isstring(s));

p.addParameter('tBounds', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1)<x(2)));
p.addParameter('smoothingFactor', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);

p.addParameter('dateLaterThan', [], @(s) isempty(s) || ischar(s) || isstring(s));
p.addParameter('dateEarlierThan', [], @(s) isempty(s) || ischar(s) || isstring(s));

% NEW
p.addParameter('day4MarkC', [], @(c) isempty(c) || iscell(c) || isstring(c));

p.addParameter('FadeToWhite', 0.85, @(x) isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('LineWidth', 2, @(x) isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MakeFigure', true, @(x) islogical(x)&&isscalar(x));
p.addParameter('lineColor', [0 0 0], @(x) isnumeric(x) && numel(x)==3 && all(x>=0 & x<=1));

p.addParameter('figSaveDir', '', @(s) isempty(s) || ischar(s) || isstring(s));

p.parse(prj_glmA, headerC, timestamps, varargin{:});
opt = p.Results;

mouseId    = string(opt.mouseId);
projType   = lower(string(opt.projType));
targetName = string(opt.targetName);

assert(strlength(mouseId)>0, 'mouseId is required (e.g., "m1045").');
assert(strlength(targetName)>0, 'targetName is required (e.g., "NoGoToneOn_1").');

% -------------------- timestamps sanity --------------------
tvec = opt.timestamps(:)'; % 1 x T

% -------------------- resolve effective dateLaterThan --------------------
% If day4MarkC provided, override dateLaterThan by mouse-specific Day4 date.
effectiveDateLaterThan = opt.dateLaterThan;

if ~isempty(opt.day4MarkC)
    d4 = opt.day4MarkC;
    if isstring(d4), d4 = cellstr(d4); end
    
    % normalize to cell array
    if iscell(d4)
        % accept either Nx2 cell OR 1D cellstr pairs - but assume Nx2 as specified
        assert(size(d4,2) == 2, 'day4MarkC must be an Nx2 cell/string array: {mouseId, "MMDDYY"; ...}.');
        mouseCol = string(d4(:,1));
        dateCol  = string(d4(:,2));
        
        hit = find(mouseCol == mouseId, 1, 'first');
        if isempty(hit)
            error('day4MarkC provided, but mouseId=%s was not found in day4MarkC.', mouseId);
        end
        
        effectiveDateLaterThan = dateCol(hit);
    else
        error('day4MarkC must be empty or an Nx2 cell/string array.');
    end
end

% -------------------- pick projection container --------------------
switch projType
    case "global"
        assert(isfield(prj_glmA,'global') && isstruct(prj_glmA.global) && isfield(prj_glmA.global,'ZC') && isfield(prj_glmA.global,'names'), ...
            'prj_glmA.global must contain ZC and names.');
        ZC_all   = prj_glmA.global.ZC;       % [J x S] cell, each [N x T x D]
        nameList = prj_glmA.global.names;    % 1 x D cellstr
    case "permouse"
        assert(isfield(prj_glmA,'perMouse') && isstruct(prj_glmA.perMouse) && isfield(prj_glmA.perMouse,'ZC') && isfield(prj_glmA.perMouse,'namesC'), ...
            'prj_glmA.perMouse must contain ZC and namesC.');
        ZC_all   = prj_glmA.perMouse.ZC;     % [J x S] cell, each [N x T x D]
        nameList = []; %#ok<NASGU>
    otherwise
        error('Unknown projType: %s (use "global" or "perMouse")', projType);
end

assert(iscell(ZC_all) && isequal(size(ZC_all), size(headerC)), 'ZC must be a cell array same size as headerC.');

% -------------------- locate row for this mouse --------------------
hdrS = cell(size(headerC));
for ii = 1:numel(headerC)
    v = headerC{ii};
    if isstring(v) && isscalar(v)
        hdrS{ii} = char(v);
    elseif ischar(v)
        hdrS{ii} = v;
    else
        hdrS{ii} = '';
    end
end
hdrS = string(hdrS);

rowHasMouse = any(contains(hdrS, mouseId, 'IgnoreCase', true), 2);
rowIdx = find(rowHasMouse, 1, 'first');
assert(~isempty(rowIdx), 'Could not find any headerC row containing mouseId=%s.', mouseId);

% perMouse names list resolution (unchanged)
if projType == "permouse"
    nameList = [];
    
    if isfield(prj_glmA,'perMouse') && isfield(prj_glmA.perMouse,'namesC') ...
            && iscell(prj_glmA.perMouse.namesC) ...
            && numel(prj_glmA.perMouse.namesC) >= rowIdx ...
            && ~isempty(prj_glmA.perMouse.namesC{rowIdx})
        nameList = prj_glmA.perMouse.namesC{rowIdx};
    end
    
    if isempty(nameList) && isfield(prj_glmA,'perMouse') && isfield(prj_glmA.perMouse,'namesC') ...
            && iscell(prj_glmA.perMouse.namesC)
        for ii = 1:numel(prj_glmA.perMouse.namesC)
            if ~isempty(prj_glmA.perMouse.namesC{ii})
                nameList = prj_glmA.perMouse.namesC{ii};
                break;
            end
        end
    end
    
    if isempty(nameList) && isfield(prj_glmA,'global') && isfield(prj_glmA.global,'names') ...
            && ~isempty(prj_glmA.global.names)
        nameList = prj_glmA.global.names;
    end
    
    assert(~isempty(nameList), 'Could not resolve axis names (perMouse.namesC empty and global.names missing).');
    nameList = cellstr(string(nameList(:)'));
end

% target column
targetCol = glmNameToColumns(nameList, cellstr(targetName));
assert(~isempty(targetCol), 'targetName=%s not found in %s names.', targetName, projType);
targetCol = targetCol(1);

% -------------------- gather sessions for that row --------------------
hdrRow = headerC(rowIdx, :);
ZC_row = ZC_all(rowIdx, :);

% optional trial selectors
useTrials = false;
if ~isempty(opt.trIdC)
    assert(iscell(opt.trIdC) && isequal(size(opt.trIdC), size(headerC)), 'trIdC must be same size as headerC.');
    trIdRow = opt.trIdC(rowIdx, :);
    useTrials = true;
else
    trIdRow = cell(size(hdrRow));
end

% -------------------- session date parsing + filtering --------------------
sessDt = NaT(1, numel(hdrRow));
sessOk = false(1, numel(hdrRow));
for jj = 1:numel(hdrRow)
    h = hdrRow{jj};
    if isempty(h), continue; end
    [dtSess, ok] = parse_header_mmddyy_ignore_suffix(h);
    if ok
        sessDt(jj) = dtSess;
        sessOk(jj) = true;
    end
end

% Later-than filter (effectiveDateLaterThan may be empty)
if ~isempty(effectiveDateLaterThan)
    dt0 = datetime(char(string(effectiveDateLaterThan)), 'InputFormat','MMddyy');
    keepLater = (sessDt >= dt0);
else
    keepLater = true(size(sessDt));
end

% Earlier-than filter (unchanged)
if ~isempty(opt.dateEarlierThan)
    dt1 = datetime(char(string(opt.dateEarlierThan)), 'InputFormat','MMddyy');
    keepEarlier = (sessDt <= dt1);
else
    keepEarlier = true(size(sessDt));
end

keepMask = ~cellfun(@isempty, hdrRow) & ~cellfun(@isempty, ZC_row) & keepLater & keepEarlier;
idxKeep = find(keepMask);
assert(~isempty(idxKeep), 'No sessions survived selection for mouseId=%s.', mouseId);

% Sort by date when parseable
dtKeep = sessDt(idxKeep);
if all(~isnat(dtKeep))
    [~, ord] = sort(dtKeep, 'ascend');
    idxKeep = idxKeep(ord);
    dtKeep = dtKeep(ord);
else
    dtKeep = sessDt(idxKeep);
end

% -------------------- time bounds --------------------
firstZ = [];
for jj = idxKeep
    if ~isempty(ZC_row{jj})
        firstZ = ZC_row{jj};
        break;
    end
end
assert(~isempty(firstZ), 'Internal: no non-empty Z found after filtering.');
assert(ndims(firstZ)==3, 'Each Z must be N x T x D.');

T = size(firstZ,2);
assert(numel(tvec)==T, 'timestamps length (%d) must match Z time dimension T (%d).', numel(tvec), T);

tMask = true(1, T);
if ~isempty(opt.tBounds)
    tMask = (tvec >= opt.tBounds(1)) & (tvec <= opt.tBounds(2));
end
tIdx = find(tMask);
assert(~isempty(tIdx), 'tBounds excluded all timestamps.');
tPlot = tvec(tIdx);

% -------------------- plot --------------------
if opt.MakeFigure
    hFig = figure('Color','w'); hold on;
else
    hFig = [];
end

Ns = numel(idxKeep);
zMeanC = cell(1, Ns);
headersUsed = cell(1, Ns);

wVec = linspace(opt.FadeToWhite, 0, Ns);
baseCol = opt.lineColor;

hLeg = gobjects(1, Ns);
legC = cell(1, Ns);

for k = 1:Ns
    jj = idxKeep(k);
    
    Z = ZC_row{jj};
    hdrSess = string(hdrRow{jj});
    headersUsed{k} = char(hdrSess);
    
    z1 = squeeze(Z(:, :, targetCol));   % [N x T]
    N  = size(z1,1);
    
    if useTrials
        trS = trIdRow{jj};
        if ~isempty(trS) && isstruct(trS) && isfield(trS, char(string(opt.trialField)))
            trI = trS.(char(string(opt.trialField)));
            if ~islogical(trI)
                tmp = false(N,1);
                tmp(trI(:)) = true;
                trI = tmp;
            end
            if numel(trI)==N && any(trI)
                z1 = z1(trI, :);
            end
        end
    end
    
    zMean = mean(z1, 1, 'omitnan');   % 1 x T
    zMean = zMean(tIdx)';             % [Tsel x 1]
    
    if opt.smoothingFactor > 0
        zMean = double(zMean);
        zMean = smooth2a(zMean, opt.smoothingFactor, 0);
    end
    zMeanC{k} = zMean;
    
    col = blend_to_white(baseCol, wVec(k));
    
    if opt.MakeFigure
        hLeg(k) = plot(tPlot, zMean, 'LineWidth', opt.LineWidth, 'Color', col);
    end
    legC{k} = char(hdrSess);
end

if opt.MakeFigure
    xlabel('time', 'Interpreter', 'none');
    ylabel(sprintf('proj (%s)', targetName), 'Interpreter', 'none');
    title(sprintf('%s | %s | %s | %s', mouseId, projType, targetName, char(string(opt.trialField))), 'Interpreter','none');
    legend(hLeg(isgraphics(hLeg)), legC(isgraphics(hLeg)), 'Location','eastoutside', 'Interpreter','none');
    box off;
end

% -------------------- save figure --------------------
figSaveDir = char(string(opt.figSaveDir));
if opt.MakeFigure && ~isempty(figSaveDir)
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end
    
    todayStr = datestr(now, 'mmddyy');
    trialFieldStr = char(string(opt.trialField));
    projTypeStr   = char(projType);
    targetStr     = char(targetName);
    
    fbase = sprintf('%s_%s_%s_%s_acrossSession_%s.pdf', ...
        char(mouseId), trialFieldStr, projTypeStr, targetStr, todayStr);
    
    fbase = sanitize_filename(fbase);
    fpath = fullfile(figSaveDir, fbase);
    
    set(gcf, 'PaperPositionMode','auto');
    print(gcf, fpath, '-dpdf', '-painters', '-bestfit');
    
    fprintf('[perMouseAcrossSession] saved: %s\n', fpath);
end

% -------------------- pack output --------------------
out = struct();
out.mouseId      = char(mouseId);
out.projType     = char(projType);
out.targetName   = char(targetName);
out.targetCol    = targetCol;
out.headersUsed  = headersUsed;
out.sessDtUsed   = dtKeep(:);
out.t            = tPlot(:);
out.zMeanC       = zMeanC;
out.figHandle    = hFig;

end

%% ===================== HELPERS (keep at end) =====================

function cols = glmNameToColumns(glmNameListC, targetNameC)
glmNameS = string(glmNameListC(:));
targS    = string(targetNameC(:));
cols = [];
for i = 1:numel(targS)
    hit = find(glmNameS == targS(i), 1, 'first');
    if ~isempty(hit)
        cols(end+1) = hit; %#ok<AGROW>
    end
end
end

function [dtSess, ok] = parse_header_mmddyy_ignore_suffix(header)
%PARSE_HEADER_MMDDYY_IGNORE_SUFFIX
% Supports "m####_MMDDYY" and "m####_MMDDYY-1", but ignores the "-1" part.
h = char(string(header));
tok = regexp(h, '_(\d{6})(?:-\d+)?$', 'tokens', 'once'); % ignore suffix
if isempty(tok)
    dtSess = NaT; ok = false; return;
end

mmddyy = tok{1};
try
    dtSess = datetime(mmddyy, 'InputFormat','MMddyy');
    ok = true;
catch
    dtSess = NaT; ok = false;
end
end

function colOut = blend_to_white(colIn, w)
colIn = colIn(:)';
if numel(colIn) ~= 3, colIn = [0 0 0]; end
w = max(min(w,1),0);
colOut = (1-w)*colIn + w*[1 1 1];
colOut = max(min(colOut,1),0);
end

function fn = sanitize_filename(fn)
fn = char(fn);
bad = '<>:"/\|?*';
for k = 1:numel(bad)
    fn(fn==bad(k)) = '_';
end
fn = regexprep(fn, '\s+', '_');
end
