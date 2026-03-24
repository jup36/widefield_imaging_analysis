function [gName, pcIdx] = parse_axis_name(nm)
nm = char(string(nm));
tok = regexp(nm, '^(.*)_(\d+)$', 'tokens', 'once');
if isempty(tok)
    gName = nm;
    pcIdx = 1;
else
    gName = tok{1};
    pcIdx = str2double(tok{2});
    if ~isfinite(pcIdx) || pcIdx < 1, pcIdx = 1; end
end
end
