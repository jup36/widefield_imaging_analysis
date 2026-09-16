function savePdf(figH, pathNoExt)
outFile = [char(pathNoExt) '.pdf'];
set(figH, 'InvertHardcopy', 'off');
print(figH, outFile, '-dpdf', '-painters');
fprintf('Saved figure:\n  %s\n', outFile);
end