function s = sanitize_for_slurm(s)
% Make a string safe for Slurm job names and script names.

s = char(s);
s = regexprep(s, '[^\w\-]', '_');

% Slurm job names can get annoying if very long
maxLen = 80;
if numel(s) > maxLen
    s = s(1:maxLen);
end

end