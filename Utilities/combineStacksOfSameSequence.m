function combineStacksOfSameSequence(imgStackC, file_list_imgStack)

assert(length(imgStackC)==length(file_list_imgStack))
folderC = cellfun(@fileparts, file_list_imgStack, 'UniformOutput', false); 
uniqueFolderC = unique(folderC); 

combinedImgStackC = cell(1, length(uniqueFolderC)); 

for i = 1:length(uniqueFolderC) 
    

end


end