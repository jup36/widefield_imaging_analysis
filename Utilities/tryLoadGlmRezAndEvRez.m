function [glmRez, glmEvRez] = tryLoadGlmRezAndEvRez(filePathBase, mId, header, glmFileKeyword)
% Locate and load one session's glmRez (for X_names) and glmEvRez (for
% cvR2_unique), mirroring the same mId -> header -> task -> Matfiles
% discovery chain used elsewhere in this pipeline. Returns [], [] (not an
% error) if anything in the chain is missing -- caller logs the reason.
glmRez = [];
glmEvRez = [];

animalFolderC = find_keyword_containing_folder(filePathBase, mId, 'recursive', false);
if isempty(animalFolderC), return; end

sessionFolderC = find_keyword_containing_folder(animalFolderC{1}, header, 'recursive', false);
if isempty(sessionFolderC), return; end

taskFolderC = find_keyword_containing_folder(sessionFolderC{1}, 'task', 'recursive', false);
if isempty(taskFolderC), return; end

filePath_mat = cell2mat(find_keyword_containing_folder(taskFolderC{1}, 'Matfiles', 'recursive', false));
if isempty(filePath_mat), return; end

glmFileC = find_keyword_containing_files(filePath_mat, glmFileKeyword, 'recursive', false);
if isempty(glmFileC), return; end

if numel(glmFileC) > 1
    fileInfo = cellfun(@(f) dir(f), glmFileC);
    [~, mostRecentIdx] = max([fileInfo.datenum]);
    glmFile = glmFileC{mostRecentIdx};
else
    glmFile = glmFileC{1};
end

try
    S = load(glmFile, 'glmRez', 'glmEvRez');
    glmRez = S.glmRez;
    glmEvRez = S.glmEvRez;
catch
    glmRez = [];
    glmEvRez = [];
end
end