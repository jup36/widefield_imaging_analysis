function mouseId = inferMouseIdFromRow_(hdrRow)
% Get mouse ID from first non-empty header in a row.

mouseId = '';
for k = 1:numel(hdrRow)
    h = hdrRow{k};
    if isempty(h), continue; end
    tok = regexp(string(h), '(m\d{3,5})', 'tokens', 'once');
    if ~isempty(tok)
        mouseId = char(string(tok{1}));
        return;
    end
end
end
