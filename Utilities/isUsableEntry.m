function tf = isUsableEntry(entry, extractFn)
tf = false;
if isempty(entry) || ~isstruct(entry), return; end
try
    M = extractFn(entry);
catch
    return;   % field missing
end
tf = ~isempty(M);
end
