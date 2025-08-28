function outLogic = verifyTrainTestSets(filePath_fnC)
load(filePath_fnC, 'fnC')
fileNames = cellfun(@getFileNameOnly, fnC, 'UniformOutput', false);
pattern = '_([0-9]{1,2})_';
numberStr = cellfun(@(f) cell2mat(regexp(f, pattern, 'tokens', 'once')), fileNames, 'UniformOutput', false);
numberVal = cellfun(@(x) str2double(x), numberStr);

% To assert first column is odd second column equals to 'first column + 1'
outLogic = all(mod(numberVal(:,1), 2) == 1) & isequal(numberVal(:,1) + 1, numberVal(:,2));

    function fname = getFileNameOnly(fpath)
        [~, name, ext] = fileparts(fpath);
        fname = [name, ext];
    end
end

