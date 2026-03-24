function [hdrC, dPMat, hdrC_all, dPMat_all] = getExpertHeadersFromdPrmC_thresholded(dPrmC, headerC, varargin)
%getExpertHeadersFromdPrmC_thresholded
% Pick up to N sessions per mouse, prioritizing threshold-passing "end" sessions,
% BUT constrained to sessions that actually exist in headerC (glmRezC universe).
% If a desired end session is missing from headerC, we "slide down" to the next
% earlier available session.
%
% OUTPUTS
%   hdrC      : (nSelectMice x numberOfSessionsMax) headers, chronological (earlier->later)
%   dPMat     : (nSelectMice x numberOfSessionsMax) d' matched to hdrC (NaN for empty)
%   hdrC_all  : (nQualMice   x numberOfSessionsMax) for ALL qualifying mice
%   dPMat_all : (nQualMice   x numberOfSessionsMax)
%
% REQUIRED INPUTS
%   dPrmC     : cell (M x 2) where {m,1}=headers (cellstr/string/char), {m,2}=d' vector
%   headerC   : cell (J x S) session headers that exist in your glmRezC/trIdC universe
%
% NAME-VALUE
%   'excludeMice'          : cellstr/string of mouse IDs to exclude
%   'selectMice'           : cellstr/string of mouse IDs to keep in hdrC output (subset)
%   'dprimeThresh'         : scalar threshold (default 1.5)
%   'numberOfSessionsMax'  : integer (default 3)
%
% RULES (per mouse)
%   1) Build the candidate list from dPrmC, but KEEP ONLY headers present in headerC.
%      (This is the key change vs. referring to dPrmC alone.)
%   2) Prefer up to N sessions with d' >= thresh, taking the LAST N (end) among those.
%   3) If fewer than N passing sessions exist (after headerC filter), fill remaining
%      slots with the LAST sessions available (even if below thresh).
%   4) Outputs are in chronological order (earlier session first), padded with [] / NaN.

% -------------------- parse --------------------
p = inputParser;
p.addParameter('excludeMice', {}, @(c) iscell(c) || isstring(c));
p.addParameter('selectMice',  {}, @(c) iscell(c) || isstring(c));
p.addParameter('dprimeThresh', 1.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('numberOfSessionsMax', 3, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.parse(varargin{:});
opt = p.Results;

excludeMouseS = string(opt.excludeMice(:));
selectMouseS  = string(opt.selectMice(:));
dprimeThresh  = double(opt.dprimeThresh);
Nmax          = double(opt.numberOfSessionsMax);

% -------------------- build availability set from headerC --------------------
hdrAvail = headerC(:);
hdrAvail = hdrAvail(~cellfun(@isempty, hdrAvail));
hdrAvailS = unique(string(hdrAvail), 'stable'); %#ok<NASGU>
availSet = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(hdrAvailS)
    availSet(char(hdrAvailS(i))) = true;
end

% -------------------- iterate dPrmC rows --------------------
all_mouseS  = strings(0,1);
all_hdrList = {};   % each entry: cellstr (1..Nmax)
all_dpList  = {};   % each entry: double  (1..Nmax)

M = size(dPrmC,1);

for m = 1:M
    if size(dPrmC,2) < 2, continue; end

    hdrCell = dPrmC{m,1};
    dp      = dPrmC{m,2};

    if isempty(hdrCell) || isempty(dp), continue; end

    % Normalize header list
    if ischar(hdrCell) || isstring(hdrCell)
        hdrCell = cellstr(hdrCell);
    end
    hdrS = string(hdrCell(:));
    dp   = double(dp(:));

    good = isfinite(dp) & (strlength(hdrS) > 0);
    hdrS = hdrS(good);
    dp   = dp(good);
    if isempty(dp), continue; end

    % Infer mouseId from first header
    tok0 = regexp(char(hdrS(1)), '(m\d{3,5})', 'tokens', 'once');
    if isempty(tok0), continue; end
    mouseId = string(tok0{1});

    if any(mouseId == excludeMouseS), continue; end

    % Restrict to headers that actually belong to this mouseId
    sameMouse = false(numel(hdrS),1);
    for ii = 1:numel(hdrS)
        tokEach = regexp(char(hdrS(ii)), '(m\d{3,5})', 'tokens', 'once');
        if ~isempty(tokEach)
            sameMouse(ii) = (string(tokEach{1}) == mouseId);
        end
    end
    hdrS = hdrS(sameMouse);
    dp   = dp(sameMouse);
    if isempty(dp), continue; end

    % IMPORTANT: restrict to sessions that exist in headerC (glmRezC universe)
    inAvail = false(numel(hdrS),1);
    for ii = 1:numel(hdrS)
        inAvail(ii) = isKey(availSet, char(hdrS(ii)));
    end
    hdrS = hdrS(inAvail);
    dp   = dp(inAvail);
    if isempty(dp)
        % nothing to return for this mouse in glm universe -> skip mouse entirely
        continue;
    end

    % Sort chronologically
    [dtKey, okParse] = parseHeaderDatetimeWithSuffix_(hdrS);
    if okParse
        [~, ord] = sort(dtKey, 'ascend');
        hdrS = hdrS(ord);
        dp   = dp(ord);
    else
        warning('Mouse %s: could not parse dates from headers; using dPrmC order.', mouseId);
    end

    % Select indices (prefer last N among passers; fill from end if needed)
    passIdx = find(dp >= dprimeThresh);
    useIdx  = [];

    if ~isempty(passIdx)
        % take last N from passers
        take = passIdx(max(1, numel(passIdx)-Nmax+1):end);
        useIdx = take(:)';
    end

    if numel(useIdx) < Nmax
        % fill remaining with last sessions (even if below thresh), skipping duplicates
        allIdxRev = numel(dp):-1:1; % from end backwards
        for ii = allIdxRev
            if numel(useIdx) >= Nmax, break; end
            if any(useIdx == ii), continue; end
            useIdx(end+1) = ii; %#ok<AGROW>
        end
    end

    % Now enforce chronological order in output (earlier first)
    useIdx = unique(useIdx, 'stable');
    useIdx = sort(useIdx, 'ascend');

    hdrOut = hdrS(useIdx);
    dpOut  = dp(useIdx);

    % If fewer than Nmax available (should be rare), pad
    if numel(hdrOut) < Nmax
        padN = Nmax - numel(hdrOut);
        hdrOut = [hdrOut; strings(padN,1)];
        dpOut  = [dpOut;  nan(padN,1)];
    elseif numel(hdrOut) > Nmax
        % just in case (should not happen), keep last Nmax but chronological
        hdrOut = hdrOut(end-Nmax+1:end);
        dpOut  = dpOut(end-Nmax+1:end);
    end

    % Convert blank strings -> [] for hdr cell output
    hdrCellOut = cell(1,Nmax);
    for k = 1:Nmax
        if strlength(hdrOut(k))==0
            hdrCellOut{k} = [];
        else
            hdrCellOut{k} = char(hdrOut(k));
        end
    end

    all_mouseS(end+1,1)  = mouseId;                 %#ok<AGROW>
    all_hdrList{end+1,1} = hdrCellOut;              %#ok<AGROW>
    all_dpList{end+1,1}  = dpOut(:)';               %#ok<AGROW>
end

% -------------------- pack ALL qualifying outputs --------------------
hdrC_all   = padHeaderRowListToCellMatrix_(all_hdrList, Nmax);
dPMat_all  = padDpRowListToNumericMatrix_(all_dpList, Nmax);

% -------------------- apply selectMice --------------------
if isempty(selectMouseS)
    keep = true(size(all_mouseS));
else
    keep = ismember(all_mouseS, selectMouseS);
end

hdrC  = hdrC_all(keep,:);
dPMat = dPMat_all(keep,:);

end

%% ========================= helpers =========================

function [dtKey, ok] = parseHeaderDatetimeWithSuffix_(hdrS)
% Parse headers like m1613_050725 or m1613_050725-1 into sortable datetimes
hdrS = string(hdrS(:));
n = numel(hdrS);
dtKey = NaT(n,1);
ok = true;

for i = 1:n
    h = char(hdrS(i));
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
    dtKey(i) = d0 + seconds(suf); % within-day ordering
end
end

function C = padHeaderRowListToCellMatrix_(hdrRowList, Nmax)
% hdrRowList: each entry is 1xNmax cell with [] or char
nMouse = numel(hdrRowList);
C = cell(nMouse, Nmax);
for i = 1:nMouse
    row = hdrRowList{i};
    if isempty(row), continue; end
    row = row(:)';
    L = min(numel(row), Nmax);
    C(i,1:L) = row(1:L);
end
end

function M = padDpRowListToNumericMatrix_(dpRowList, Nmax)
nMouse = numel(dpRowList);
M = nan(nMouse, Nmax);
for i = 1:nMouse
    row = dpRowList{i};
    if isempty(row), continue; end
    row = double(row(:)');
    L = min(numel(row), Nmax);
    M(i,1:L) = row(1:L);
end
end