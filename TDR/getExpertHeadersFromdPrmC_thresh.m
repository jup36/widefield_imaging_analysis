function [hdrC, dpBest, hdrC_all] = getExpertHeadersFromdPrmC_thresh(dPrmC, varargin)
%GETEXPERTHEADERSFROMDPRMC_THRESH
% Return best-d' session header per mouse for mice with at least one session
% meeting a d' threshold. Optionally (multiSession=true) also return later
% sessions for that mouse (only if later sessions also pass the threshold).
%
% OUTPUTS
%   hdrC      : cell array (nSelectMice x nSessMax) AFTER applying selectMice (if provided)
%               Each row = one mouse, each column = one header string (or empty []).
%   dpBest    : numeric vector (nSelectMice x 1), best d' per selected mouse (among passing)
%   hdrC_all  : cell array (nQualMice x nSessMaxAll), best-per-mouse (+ later sessions if multiSession)
%               for ALL qualifying mice (exclude + thresh only)
%
% NAME-VALUE
%   'excludeMice'   : cellstr/string of mouse IDs to exclude (e.g., {'m1893'})
%   'selectMice'    : cellstr/string of mouse IDs to keep in hdrC output (subset)
%   'dprimeThresh'  : scalar threshold
%   'multiSession'  : false (default). If true, include later sessions after best
%                     (within that mouse only) that also pass threshold.
%
% NOTES
% - Headers are assumed to contain mouseId like m#### or m##### and date like _MMDDYY.
% - Supports same-day repeats like: m1613_050725-1 (sorted after m1613_050725).
% - Robust against the bug where sessions from other mice leak into a mouse row:
%   we explicitly filter headers by mouseId before any multiSession selection.

% -------------------- parse --------------------
p = inputParser;
p.addParameter('excludeMice', {}, @(c) iscell(c) || isstring(c));
p.addParameter('selectMice',  {}, @(c) iscell(c) || isstring(c));
p.addParameter('dprimeThresh', -Inf, @(x) isnumeric(x) && isscalar(x));
p.addParameter('multiSession', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

excludeMouseS = string(opt.excludeMice(:));
selectMouseS  = string(opt.selectMice(:));
dprimeThresh  = opt.dprimeThresh;
multiSession  = opt.multiSession;

% -------------------- collect per-mouse lists (ALL qualifying) --------------------
all_mouseS   = strings(0,1);     % (nMouse x 1)
all_hdrList  = {};              % cell, each entry is cellstr list for that mouse (1 x nSess)
all_dpBest   = [];              % best d' per mouse (scalar)

M = size(dPrmC,1);

for m = 1:M
    if size(dPrmC,2) < 2
        continue;
    end

    hdrCell = dPrmC{m,1};
    dp      = dPrmC{m,2};

    if isempty(hdrCell) || isempty(dp)
        continue;
    end

    % Normalize hdrCell -> cellstr
    if ischar(hdrCell) || isstring(hdrCell)
        hdrCell = cellstr(hdrCell);
    end

    hdrS = string(hdrCell(:));
    dp   = double(dp(:));

    % Valid entries
    good = isfinite(dp) & (strlength(hdrS) > 0);
    hdrS = hdrS(good);
    dp   = dp(good);
    if isempty(dp)
        continue;
    end

    % Infer mouseId from first header
    tok0 = regexp(char(hdrS(1)), '(m\d{3,5})', 'tokens', 'once');
    if isempty(tok0)
        continue;
    end
    mouseId = string(tok0{1});

    % Exclude
    if any(mouseId == excludeMouseS)
        continue;
    end

    % ---- CRITICAL: restrict to headers that actually match this mouseId ----
    sameMouse = false(numel(hdrS),1);
    for ii = 1:numel(hdrS)
        tokEach = regexp(char(hdrS(ii)), '(m\d{3,5})', 'tokens', 'once');
        if ~isempty(tokEach)
            sameMouse(ii) = (string(tokEach{1}) == mouseId);
        end
    end
    hdrS = hdrS(sameMouse);
    dp   = dp(sameMouse);
    if isempty(dp)
        continue;
    end

    % Threshold qualification
    pass = dp >= dprimeThresh;
    if ~any(pass)
        continue;
    end

    % Sort sessions for THIS mouse: MMDDYY, then optional "-k" suffix
    [dtKey, okParse] = parseHeaderDatetimeWithSuffix_(hdrS);

    if okParse
        [~, ord] = sort(dtKey, 'ascend');
        hdrS = hdrS(ord);
        dp   = dp(ord);
        pass = pass(ord);
    else
        warning('Mouse %s: could not parse MMDDYY from headers; multiSession ordering fell back to original header order.', mouseId);
    end

    % Best among passing (for this mouse)
    dpPass  = dp(pass);
    hdrPass = hdrS(pass);
    [bestDp, imx] = max(dpPass);
    bestHdr = hdrPass(imx);

    % Build output header list for this mouse
    if ~multiSession
        outHdrS = bestHdr;
    else
        % Find the best header position in the FULL (sorted) list (not just pass list)
        iBest = find(hdrS == bestHdr, 1, 'first');

        outHdrS = bestHdr;
        if ~isempty(iBest)
            laterIdx = ((1:numel(hdrS))' > iBest) & pass;  % later AND passes threshold
            outHdrS = [outHdrS; hdrS(laterIdx)];
        end

        % De-duplicate (just in case)
        outHdrS = unique(outHdrS, 'stable');
    end

    % Store
    all_mouseS(end+1,1)  = mouseId;            %#ok<AGROW>
    all_hdrList{end+1,1} = cellstr(outHdrS);   %#ok<AGROW>
    all_dpBest(end+1,1)  = bestDp;             %#ok<AGROW>
end

% -------------------- pack hdrC_all as padded cell array --------------------
hdrC_all = padHeaderListToCellMatrix_(all_hdrList);

% -------------------- apply selectMice to form hdrC / dpBest --------------------
if isempty(selectMouseS)
    keep = true(size(all_mouseS));
else
    keep = ismember(all_mouseS, selectMouseS);
end

hdrC   = padHeaderListToCellMatrix_(all_hdrList(keep));
dpBest = all_dpBest(keep);

end

%% ========================= helpers (end of file) =========================
function [dtKey, ok] = parseHeaderDatetimeWithSuffix_(hdrS)
% Parse headers like:
%   m1613_050725
%   m1613_050725-1
% into sortable datetime keys.
%
% Rule: same MMDDYY, suffix -k sorts AFTER -0 (k=0 for no suffix).

hdrS = string(hdrS(:));
n = numel(hdrS);
dtKey = NaT(n,1);
ok = true;

for i = 1:n
    h = char(hdrS(i));
    % capture _MMDDYY or _MMDDYY-k at end
    tok = regexp(h, '_(\d{6})(?:-(\d+))?$', 'tokens', 'once');
    if isempty(tok)
        ok = false;
        return;
    end

    mmddyy = tok{1};
    if numel(tok) >= 2 && ~isempty(tok{2})
        suf = str2double(tok{2});
        if ~isfinite(suf), suf = 0; end
    else
        suf = 0;
    end

    d0 = datetime(mmddyy, 'InputFormat','MMddyy');
    dtKey(i) = d0 + seconds(suf); % tiny offset for within-day ordering
end
end

function C = padHeaderListToCellMatrix_(hdrList)
% hdrList: cell array where each element is a cellstr vector (sessions for one mouse)
% returns padded cell matrix (nMouse x nSessMax), empties are [].

if isempty(hdrList)
    C = cell(0,0);
    return;
end

nMouse = numel(hdrList);
lens = zeros(nMouse,1);
for i = 1:nMouse
    if isempty(hdrList{i})
        lens(i) = 0;
    else
        lens(i) = numel(hdrList{i});
    end
end
nSessMax = max(lens);

C = cell(nMouse, nSessMax);
for i = 1:nMouse
    if lens(i) == 0, continue; end
    C(i,1:lens(i)) = hdrList{i}(:)'; % row fill
end
end
