function [expertHeaderC, expertMouseC, expertDprime] = getExpertHeadersFromdPrmC(dPrmC, excludeMouseC)
%GETEXPERTHEADERSFROMDPRMC  Return best-d' session header per mouse.
%
% INPUT
%   dPrmC : Mx3 cell array
%           dPrmC{m,1} = 1xS cell array of headers (e.g., 'm1045_122424')
%           dPrmC{m,2} = 1xS numeric array of d' values
%   excludeMouseC : cellstr of mouse IDs to ignore (e.g., {'m1893',...})
%
% OUTPUT
%   expertHeaderC : {Nmice x 1} best session header per mouse
%   expertMouseC  : {Nmice x 1} mouse IDs
%   expertDprime  : [Nmice x 1] best d' per mouse

if nargin < 2 || isempty(excludeMouseC)
    excludeMouseC = {};
end
excludeMouseC = string(excludeMouseC(:));

M = size(dPrmC,1);

expertHeaderC = {};
expertMouseC  = {};
expertDprime  = [];

for m = 1:M
    if size(dPrmC,2) < 2, continue; end
    hdrC = dPrmC{m,1};
    dp   = dPrmC{m,2};

    if isempty(hdrC) || isempty(dp), continue; end

    % Normalize header cell
    if isstring(hdrC) || ischar(hdrC)
        hdrC = cellstr(hdrC);
    end
    hdrS = string(hdrC(:));
    dp   = double(dp(:));

    % Drop invalid
    good = isfinite(dp) & (strlength(hdrS) > 0);
    hdrS = hdrS(good);
    dp   = dp(good);
    if isempty(dp), continue; end

    % Infer mouse ID from first header
    tok = regexp(hdrS(1), '(m\d{3,5})', 'tokens', 'once');
    if isempty(tok)
        % fallback: try all headers
        tok = regexp(hdrS, '(m\d{3,5})', 'tokens', 'once');
        tok = tok(~cellfun(@isempty,tok));
        if isempty(tok), continue; end
        mouseId = string(tok{1}{1});
    else
        mouseId = string(tok{1});
    end

    % Exclude mice if requested
    if any(mouseId == excludeMouseC)
        continue;
    end

    % Best d' (ties -> first)
    [mx, imx] = max(dp);
    bestHdr = hdrS(imx);

    expertMouseC{end+1,1}  = char(mouseId);
    expertHeaderC{end+1,1} = char(bestHdr);
    expertDprime(end+1,1)  = mx; %#ok<AGROW>
end
end
