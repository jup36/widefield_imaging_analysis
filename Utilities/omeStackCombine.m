function omeStackCombine(filePath)

% directory for imaging data
fileImg = GrabFiles_sort_trials('_img', 0, {fullfile(filePath)});

% load preprocessed stacks ('ome_stack.mat')
[file_list_img, ~] = GrabFiles_sort_trials('ome_stack.mat', 1, fileImg(1)); % use GrabFiles_sort_trials to sort both files and folders

folderC = cellfun(@fileparts, file_list_img, 'UniformOutput', false); % folder cell (the ome_stack files in the same folder will be combined)
uniqueFolderC = unique(folderC); 
%Reorder by the natural order of the file names (e.g. 1,2,3 not 1,10,100)
[~, uniqueFolderI]=natsortfiles(uniqueFolderC);
uniqueFolderC  = uniqueFolderC(uniqueFolderI);

% combine and save ome stacks from the same folder (sequence)
for i = 1:length(uniqueFolderC) 
    filesIdxForThisFolder = find(cellfun(@(a) strcmpi(a, uniqueFolderC{i}), folderC)); 
    imgC = cell(1, length(filesIdxForThisFolder)); 
   
    for ii = 1:length(filesIdxForThisFolder)
        load(file_list_img{1, filesIdxForThisFolder(ii)}, 'stack')
        if ii == 1
            load(file_list_img{1, filesIdxForThisFolder(ii)}, 'opts'); % load stacks to be combined
        end
        imgC{1, ii} = stack; 
        fprintf("finished loading stack file #%d\n", filesIdxForThisFolder(ii))
    end

    % sanity check to ensure that the image stacks are ordered properly and have equal dimensions 
    assert(length(unique(cellfun(@(a) size(a, 1), imgC)))==1) % the number of rows match
    assert(length(unique(cellfun(@(a) size(a, 2), imgC)))==1) % the number of columns match
    assert(~any(diff(cellfun(@(a) size(a, 3), imgC))>0)) % the number of frames must descend 

    combinedStackC = cell(1, 1, length(imgC)); 
    for ii = 1:length(imgC)
        combinedStackC{:,:,ii} = imgC{1, ii}; 
    end 
    
    combinedStack = cell2mat(combinedStackC); 

    % verify BV alternation and create bvFrameI, in which 1 indicates blue
    % (brighter) frame
    bvFrameI = zeros(1, size(combinedStack, 3)); 
    mCombinedStack = squeeze(nanmean(nanmean(combinedStack, 1), 2)); 
    diffMeanLogic = diff(mCombinedStack)>0; 
    alterLogic = all(abs(diff(diffMeanLogic))==1); % if there's a frame whose luminosity doesn't alternate, there would be 0. 
    if alterLogic % frames are proven to be alternating
        if diffMeanLogic(1)==false % the combinedStack starts with brighter (B) frame
            bvFrameI(1, 1:2:end) = 1; 
        elseif diffMeanLogic(1)==true % the combinedStack starts with darker (V) frame
            bvFrameI(1, 2:2:end) = 1; 
        end
    else
        warning("There are frames not alternating correctly!") 
        bvFrameI = nan(1, size(combinedStack, 3)); ss
    end
    
    bvFrameIC{1, i} = bvFrameI; 

    % save
    [~, fileName] = fileparts(file_list_img{filesIdxForThisFolder(1)});
    matches = regexp(fileName, '^(.*?MMStack_)', 'tokens');
    saveNameHeader = matches{1}{1};  % The captured string includes 'MMStack_'
    saveName = [saveNameHeader, 'combinedStack']; 
    save(fullfile(uniqueFolderC{i}, saveName), 'combinedStack', 'opts', 'bvFrameI')
        
    clearvars imgC combinedStack combinedStackC bvFrameI
    fprintf("Combined and saved ome stacks in folder #%d\n", i); 
end


end