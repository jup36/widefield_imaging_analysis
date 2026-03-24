function hdrMat = normalizeHeaderMatrix_(hdrIn)
% Convert input (cell or string array) to cell matrix of char headers; empty -> [].

if isempty(hdrIn)
    hdrMat = cell(0,0);
    return;
end

if isstring(hdrIn)
    hdrIn = cellstr(hdrIn);
end

% If user passed a vector, make it a column.
if iscell(hdrIn) && isvector(hdrIn)
    hdrIn = hdrIn(:);
end

hdrMat = cell(size(hdrIn));
for i = 1:numel(hdrIn)
    x = hdrIn{i};

    % common "empty cell shows as 0x0 double" case
    if isempty(x)
        hdrMat{i} = [];
        continue;
    end

    if isstring(x) || ischar(x)
        s = char(string(x));
        if strlength(string(s))==0
            hdrMat{i} = [];
        else
            hdrMat{i} = s;
        end
    else
        % anything else (e.g. numeric empty) -> treat as empty
        hdrMat{i} = [];
    end
end
end