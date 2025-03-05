function [data_norm, nanpxs, data_train, data_test] = ProcessAndSplitDataAuditoryGng(fnC_fn, save_fn, parameter_class)
%Camden MacDowell - timeless (Modified by Junchol Park Feb 2025)
%filters, normalizes, and splits data in fn into training and test test
%fn can be the full file path or a stack and opt structure. 
% Note for modification (Feb 2025 Junchol Park)
% This modification is to accommodate the data acquisition scheme for
%   auditory go/no-go task. 
% Input: 
%  fnC_fn: A key input to load fnC
%  fnC: A key cell array that contains directories for combined dff matfiles.
%   Importantly fnC must have two columns with the first column specifies
%   files that comprise train splits, whereas the second column specifies
%   files thta comprise test splits for fpCNMF. Thus, the num_chunks 
%   must equal to size(fnC, 1), and size(fnC, 2) must be 2. 
% Output: 
%   data_norm
%   nanpxs
%   data_train
%   data_test

load(fnC_fn, 'fnC'); 

if ~ispc
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'))
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF/'));
end

gp = loadobj(feval(parameter_class)); 

num_chunks = size(fnC, 1); 
assert(size(fnC, 2)==2); % train and test
assert(sum(cellfun(@isempty, fnC(:)))==0); 

%% Load dff stacks
for rr = 1:size(fnC, 1)
    for cc = 1:size(fnC, 2)
        temp = load(fnC{rr, cc});
        if ~exist('opts', 'var')
            opts = temp.opts; 
        end      
        % dimension
        xC{rr, cc} = size(temp.dff, 1); 
        yC{rr, cc} = size(temp.dff, 2);  
        zC{rr, cc} = size(temp.dff, 3);  
        % condition data and remove nan pxls
        [dataC{rr, cc}, nanpxsC{rr, cc}] = conditionDffMat(temp.dff); % dataC entries are frame-by-pixels(nonNaN)
        % deconvolution filter data (use lucric)
        fprintf('\n\tPerforming a Lucy-Goosey Deconvolution (Lucy-Richardson)\n')
        for px = 1:size(dataC{rr, cc},2)
            dataC{rr, cc}(:,px) = lucric(dataC{rr, cc}(:,px), gp.d_gamma, gp.d_smooth, gp.d_kernel);
        end
        
        clear temp; 
        fprintf(sprintf("Completed loading dff or row#%d and col#%d\n", rr, cc)); 
    end
end

%% dimension sanity check
minFrN = min(cell2mat(zC(:)));  % min number of frames
rowN = unique(cell2mat(xC(:))); % # rows: 64
assert(isscalar(rowN)); % match # rows
colN = unique(cell2mat(yC(:))); % # columns: 64
assert(isscalar(colN)); % match # columns
assert(isscalar(unique(cellfun(@nansum, nanpxsC(:))))); % match NaN pixels 
nanpxs = nanpxsC{1}; 

%% curtail frames if needed and reshape 
% curtail frames to match the number of frames in each chunk
dataC = cellfun(@(a) a(1:minFrN, :), dataC, 'UniformOutput', false); 

% back to stack orientation and stack
dataC = cellfun(@(a, b) conditionDffMat(a, b, [], [rowN, colN, minFrN]), dataC, nanpxsC, 'UniformOutput', false); 
data = cat(3, dataC{:}); % stack up all image stacks (e.g., 64 x 64 x N total frames)

%% Denoise with PCA (removed banded pixels)
if gp.w_pca_denoise
    data = DenoisePCA(data); % data: 64 x 64 x N total frames
end
[data,~] = conditionDffMat(data); % data: N total frames x (non-NaN pixels)

%% normalize to 0 to 1 
fprintf('\n\tPerforming %s normalization to %d value', gp.w_normalization_method, gp.w_norm_val);
switch gp.w_normalization_method
    case 'pixelwise' %each between x and xth pixel intensity
        data_norm = NaN(size(data));
        for px = 1:size(data,2)
            data_norm(:,px) = (data(:,px))/(prctile(data(:,px),gp.w_norm_val));
        end             
    case 'full' %normalize using the percentile of the maximum         
        data_norm = data/prctile(data(data>eps),gp.w_norm_val);          
    case 'bounded'
        data_norm = (data)/(gp.w_norm_val(2)); %normalize between zero and the upper bound     
    case 'none'
        data_norm = data;
    otherwise
        error('Unknown normalization method. Check general params')
end

%% transpose (fpCNMF operates rowwise)
data_norm = data_norm'; % non-NaN pixels X N total frames

%% create data_train, data_test
assert(isequal(size(data_norm, 2), minFrN*numel(dataC))); 

data_train = NaN(size(data_norm, 1), minFrN, size(dataC, 1)); % non-NaN pixels X frames per chunk X # of train chunks   
data_test = NaN(size(data_norm, 1), minFrN, size(dataC, 1)); % non-NaN pixels X frames per chunk X # of test chunks

count_trainsets = 0; 
count_testsets = 0; 
for i = 1:numel(dataC)
    tempDat = data_norm(:, (i-1)*minFrN+1:i*minFrN);
    if mod(i,2)==1
        count_trainsets = count_trainsets + 1; 
        data_train(:, :, count_trainsets) = tempDat; % odd chunks 
    else
        count_testsets = count_testsets + 1;
        data_test(:, :, count_testsets) = tempDat; % even chunks
    end
end

%% save off the data in the scratch directory and the nanpxs
if ~isempty(save_fn)
    fprintf('\n\tSaving data')
    save(save_fn,'data_norm','data_test','data_train','nanpxs','opts','gp','num_chunks','-v7.3')
    fprintf('\n\tDONE')
end

end















