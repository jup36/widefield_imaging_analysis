filePath_base = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed'; 
fileSaveDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData'; 
if ispc
    filePath_base = 'Z:\Rodent Data\Wide Field Microscopy\ExampleData\Preprocessed'; 
    fileSaveDir = 'Z:\Rodent Data\dualImaging_parkj\collectData'; 
end

mlist = {'m1873', 'm1045', 'm1859', 'm1613', 'm1044', 'm1048', 'm1049'}; 

parameter_class = 'general_params_dual';
gp = loadobj(feval(parameter_class));

%% Initialize outputs
wC = {}; 
nanpxsC = {}; 
idC = {}; 
skippedFiles = {};  % <–– NEW: to track failed file IDs

for f = 1:length(mlist)
    sessions = GrabFiles_sort_trials([mlist{f} '*motif'], 0, {filePath_base});  
    for ff = 1:length(sessions)
        header = extract_date_animalID_header(sessions{ff}); 
        searchStr = [header, '*green*chunk*.mat']; 
        filePath_chunk = findFilesWithPattern(sessions{ff}, searchStr); 
        for j = 1:length(filePath_chunk)
            try
                clearvars w w_train nanpxs
                load(filePath_chunk{j}, 'w', 'w_train', 'nanpxs')

                nanpxsC = [nanpxsC, {nanpxs}]; 

                if exist('w_train', 'var') && (size(w_train, 2) >= size(w, 2))
                    wC = [wC, {w_train}]; 
                else 
                    wC = [wC, {w}]; 
                end 

                [~, id] = fileparts(filePath_chunk{j});
                idC = [idC, {id}]; 
                fprintf("Finished loading chunk%d-file%d of m%d!\n", j, ff, f); 

            catch ME
                [failedDIR, failedID] = fileparts(filePath_chunk{j});
                skippedFiles = [skippedFiles, {fullfile(failedDIR, failedID)}];  % collect ID
                warning("Skipping corrupted or unreadable file:\n%s\nReason: %s", ...
                        filePath_chunk{j}, ME.message);
                continue
            end
        end
    end
end

% save
save(fullfile(fileSaveDir, 'clusterW_input_DAmotifs_072125'), 'wC', 'nanpxsC', 'idC', 'skippedFiles', '-v7.3'); % Also backup at '/Volumes/buschman/Rodent Data/dualImaging_parkj/motif_data_backup'

%% Retry error files (replace files with error from the backup ('/Volumes/buschman/Rodent Data/dualImaging_parkj/motif_data_backup'), then retry)
% stillMissing = {}; 
% for jj = 1:length(skippedFiles)
%     try
%         clearvars w w_train nanpxs
%         load(skippedFiles{jj}, 'w', 'w_train', 'nanpxs')
% 
%         nanpxsC = [nanpxsC, {nanpxs}];
% 
%         if exist('w_train', 'var') && (size(w_train, 2) >= size(w, 2))
%             wC = [wC, {w_train}];
%         else
%             wC = [wC, {w}];
%         end
% 
%         [~, id] = fileparts(skippedFiles{jj});
%         idC = [idC, {id}];
%         fprintf("Finished loading chunk%d-file%d of m%d!\n", j, ff, f);
% 
%     catch ME
%         [failedDIR, failedID] = fileparts(skippedFiles{jj});
%         stillMissing = [stillMissing, {fullfile(failedDIR, failedID)}];  % collect ID
%         warning("Skipping corrupted or unreadable file:\n%s\nReason: %s", ...
%             skippedFiles{jj}, ME.message);
%         continue
%     end
% end

%% run ClusterW (multi clustering levels)
% Optional load
% load(fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData', 'clusterW_input_CAmotifs'), 'wC', 'nanpxsC', 'idC', 'skippedFiles');

% Select motifs with the same NaN-pixel profile
% Convert each array in nanpxsC to a string representation
asStr = cellfun(@(x) mat2str(x), nanpxsC, 'UniformOutput', false);
% Use unique to find unique string representations and their counts
[uniqueStrs, ~, ic] = unique(asStr);
counts = accumarray(ic, 1);
% Find the index of the most common string
[~, mostCommonIdx] = max(counts);
% Retrieve the most common nanpxs array
mostCommonNanpxs = nanpxsC{find(strcmp(asStr, uniqueStrs{mostCommonIdx}), 1)};

is_equal_nanpxsC = cellfun(@(x) isequal(x, mostCommonNanpxs), nanpxsC);
wCV = wC(is_equal_nanpxsC); 
nanpxsCV = nanpxsC(is_equal_nanpxsC); 
save(fullfile(fileSaveDir, 'clusterW_input_DAmotifs_072125'), 'wCV', 'nanpxsCV', 'is_equal_nanpxsC', '-append'); % Also backup at '/Volumes/buschman/Rodent Data/dualImaging_parkj/motif_data_backup'

%[W_basis, kval, ovr_q, cluster_idx, cluster_idx_full, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs, shift] = ClusterW(wCV, gp, nanpxsCV{1}); 
gp.clust_removepad=1; % to remove temporal padding (10 frames before and after)
[W_basisC, W_aligned, kval, ovr_q, cluster_idxC, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs, shiftC] = ClusterW_multiLevels(wCV, gp, nanpxsCV{1}); 

% temporarily save ClusterW output
save(fullfile(fileSaveDir, 'clusterW_output_DAmotifs_072125'), 'W_basisC', 'W_aligned', 'kval', 'ovr_q', 'cluster_idxC', 'idx_knn', 'tcorr_mat', 'lag_mat', 'lags', 'nanpxs', 'shiftC', '-v7.3');

%PlotMotifs
%RefitBasisMotifs

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileList = findFilesWithPattern(baseDir, searchPattern)
%FINDFILESWITHPATTERN  Return full paths of all files that match a pattern.
%
%  fileList = findFilesWithPattern(baseDir, searchPattern)
%
%  INPUTS
%    baseDir       – directory to search (no recursion)
%    searchPattern – filename pattern with * wildcards, e.g.
%                    'm1873_041025*red*chunk*.mat'
%
%  OUTPUT
%    fileList      – cell array of full paths (char) matching the pattern.
%                    Returns {} if no file matches.
%
%  EXAMPLE
%    p = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/...
%         m1873_041025_task_day0_img_motif';
%    files = findFilesWithPattern(p,'m1873_041025*red*chunk*.mat');
%
%  NOTES
%    • Uses MATLAB's DIR, so all standard wildcard rules apply.
%    • No recursive descent; modify with **dir(fullfile(baseDir,'**',pattern))**
%      if you ever need that.
% ------------------------------------------------------------------------

arguments
    baseDir (1,1) string  {mustBeFolder}
    searchPattern (1,1) string
end

% Grab matches
info = dir(fullfile(baseDir, searchPattern));

% Build full paths
if isempty(info)
    fileList = {};
else
    fileList = fullfile({info.folder}, {info.name});
end
end
