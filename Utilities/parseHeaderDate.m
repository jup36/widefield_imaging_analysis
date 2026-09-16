
function dt = parseHeaderDate(headerCol)
% PARSEHEADERDATE
%   Parses headers like 'm1613_050725' or 'm1613_050725-1' into sortable
%   datetimes. Same format/regex convention as
%   getExpertHeadersFromdPrmC_thresholded's parseHeaderDatetimeWithSuffix_
%   (MMddyy date, optional '-N' suffix for same-day reruns, used only for
%   within-day tie-breaking ordering). headerCol: cellstr or string array.
%   Returns NaT for any header that doesn't match the expected pattern
%   (caller is responsible for excluding those from downstream ordering).
headerCol = string(headerCol(:));
n = numel(headerCol);
dt = NaT(n, 1);
for i = 1:n
    h = char(headerCol(i));
    tok = regexp(h, '_(\d{6})(?:-(\d+))?$', 'tokens', 'once');
    if isempty(tok)
        continue;   % left as NaT
    end
    mmddyy = tok{1};
    if numel(tok) >= 2 && ~isempty(tok{2})
        suf = str2double(tok{2});
        if ~isfinite(suf), suf = 0; end
    else
        suf = 0;
    end
    try
        d0 = datetime(mmddyy, 'InputFormat', 'MMddyy');
        dt(i) = d0 + seconds(suf);
    catch
        % leave as NaT
    end
end
end
