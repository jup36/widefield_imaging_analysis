function fname = getFileNameOnly(fpath)
[~, name, ext] = fileparts(fpath);
fname = [name, ext];
end
