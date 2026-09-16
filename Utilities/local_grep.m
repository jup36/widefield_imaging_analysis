function hits = local_grep(rootDir, pat)
f = dir(fullfile(rootDir, '**', '*.m'));
hits = {};
for i = 1:numel(f)
    p = fullfile(f(i).folder, f(i).name);
    txt = fileread(p);
    lines = strsplit(txt, newline);
    idx = find(contains(lines, pat));
    if isempty(idx), continue; end
    for k = idx(:)'
        fprintf('%s:%d: %s\n', p, k, strtrim(lines{k}));
    end
    hits{end+1} = p; %#ok<AGROW>
end
fprintf('\n%d file(s) matched.\n', numel(hits));
end