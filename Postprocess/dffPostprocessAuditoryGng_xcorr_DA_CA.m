function dffPostprocessAuditoryGng_xcorr_DA_CA(filePath)
%filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_121324/task';

%% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

fileBeh_green = GrabFiles_sort_trials('_green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
tbytDat_green = load(fullfile(fileBeh_green{1}), 'tbytDat'); 
tbytDat_green = tbytDat_green.('tbytDat'); 

fileBeh_red = GrabFiles_sort_trials('_red_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
tbytDat_red = load(fullfile(fileBeh_red{1}), 'tbytDat'); 
tbytDat_red = tbytDat_red.('tbytDat'); 

% load the preprocessed dffs collect
dffsmCell_green = load(fullfile(filePath, 'Matfiles', strcat(header, "_green_dff_smCollect.mat")), 'dffsmCell');
dffsmCell_green = dffsmCell_green.('dffsmCell'); 

dffsmCell_red = load(fullfile(filePath, 'Matfiles', strcat(header, '_red_dff_smCollect.mat')), 'dffsmCell');
dffsmCell_red = dffsmCell_red.('dffsmCell');

% load Allen dorsalMap
load(fullfile('/Users/jp3025/Documents/codes/WidefieldAnalysis_Musall/allenDorsalMap.mat'), 'dorsalMaps', 'motorMask', 'smotorMask', 'ssMask', 'vMask', 'rsMask');

% load Allen transformation object (the alignment parameters to Allen CCF)
folder_list_imgTrial = GrabFiles_sort_trials('_img', 0, {filePath}); 
if strcmpi(channel, 'green')
    path_transParamsAllen = find_keyword_containing_files(folder_list_imgTrial{1}, 'transParamsAllen_green', 'recursive', true); 
elseif strcmpi(channel, 'red')
    path_transParamsAllen = find_keyword_containing_files(folder_list_imgTrial{1}, 'transParamsAllen_red', 'recursive', true); 
end
load(fullfile(path_transParamsAllen{1}), 'transParams')

% file directory for trials
filePathTrials = fullfile(filePath, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

for t = 1:length(dffsmCell_green)
    %% green
    % align image stack to the AllenCCF
    dffG = alignStackToAllenKabsch(dffsmCell_green{t}, dorsalMaps.dorsalMap, transParams.tformObj); % row x col x frame
    frT_G = tbytDat_green(t).frameTrel; % store timestamps

    % apply region masks
    dffM1_G = apply2DMaskTo3DStack(dffG, motorMask);
    dffM2_G = apply2DMaskTo3DStack(dffG, smotorMask);
    dffSs_G = apply2DMaskTo3DStack(dffG, ssMask);
    dffV1_G = apply2DMaskTo3DStack(dffG, vMask);
    dffRs_G = apply2DMaskTo3DStack(dffG, rsMask);

    % Align to tone Onset
    [rezG.stimOnDffC.m1{t, 1}, rezG.stimOnDffC.m1{t, 2}] = alignToEvent(dffM1_G, 0, frT_G, [-0.9 5]);
    [rezG.stimOnDffC.m2{t, 1}, rezG.stimOnDffC.m2{t, 2}] = alignToEvent(dffM2_G, 0, frT_G, [-0.9 5]);
    [rezG.stimOnDffC.ss{t, 1}, rezG.stimOnDffC.ss{t, 2}] = alignToEvent(dffSs_G, 0, frT_G, [-0.9 5]);
    [rezG.stimOnDffC.v1{t, 1}, rezG.stimOnDffC.v1{t, 2}] = alignToEvent(dffV1_G, 0, frT_G, [-0.9 5]);
    [rezG.stimOnDffC.rs{t, 1}, rezG.stimOnDffC.rs{t, 2}] = alignToEvent(dffRs_G, 0, frT_G, [-0.9 5]);

    %% red
    % align image stack to the AllenCCF
    dffR = alignStackToAllenKabsch(dffsmCell_red{t}, dorsalMaps.dorsalMap, transParams.tformObj); % row x col x frame
    frT_R = tbytDat_red(t).frameTrel; % store timestamps

    % apply region masks
    dffM1_R = apply2DMaskTo3DStack(dffR, motorMask);
    dffM2_R = apply2DMaskTo3DStack(dffR, smotorMask);
    dffSs_R = apply2DMaskTo3DStack(dffR, ssMask);
    dffV1_R = apply2DMaskTo3DStack(dffR, vMask);
    dffRs_R = apply2DMaskTo3DStack(dffR, rsMask);

    % Align to tone Onset
    [rezR.stimOnDffC.m1{t, 1}, rezR.stimOnDffC.m1{t, 2}] = alignToEvent(dffM1_R, 0, frT_R, [-0.9 5]);
    [rezR.stimOnDffC.m2{t, 1}, rezR.stimOnDffC.m2{t, 2}] = alignToEvent(dffM2_R, 0, frT_R, [-0.9 5]);
    [rezR.stimOnDffC.ss{t, 1}, rezR.stimOnDffC.ss{t, 2}] = alignToEvent(dffSs_R, 0, frT_R, [-0.9 5]);
    [rezR.stimOnDffC.v1{t, 1}, rezR.stimOnDffC.v1{t, 2}] = alignToEvent(dffV1_R, 0, frT_R, [-0.9 5]);
    [rezR.stimOnDffC.rs{t, 1}, rezR.stimOnDffC.rs{t, 2}] = alignToEvent(dffRs_R, 0, frT_R, [-0.9 5]);
    
    %% xcorr
    %m1
    [rezG.dffsOnTfItp.m1{t, 1}, rezR.dffsOnTfItp.m1{t, 1}] = alignTwoDffTimeSeries(rezG.stimOnDffC.m1{t, 1}, rezR.stimOnDffC.m1{t, 1}, ...
        rezG.stimOnDffC.m1{t, 2}, rezR.stimOnDffC.m1{t, 2}, [-0.9 5]);    
    [rez.xcorrDffGR.m1{t, 1}, rez.xcorrDffGR.m1{t, 2}] = xcorr(rezG.dffsOnTfItp.m1{t, 1}-min(rezG.dffsOnTfItp.m1{t, 1}),...
        rezR.dffsOnTfItp.m1{t, 1}-min(rezR.dffsOnTfItp.m1{t, 1}), 2000, 'normalized');

    %m2
    [rezG.dffsOnTfItp.m2{t, 1}, rezR.dffsOnTfItp.m2{t, 1}] = alignTwoDffTimeSeries(rezG.stimOnDffC.m2{t, 1}, rezR.stimOnDffC.m2{t, 1}, ...
        rezG.stimOnDffC.m2{t, 2}, rezR.stimOnDffC.m2{t, 2}, [-0.9 5]);
    [rez.xcorrDffGR.m2{t, 1}, rez.xcorrDffGR.m2{t, 2}] = xcorr(rezG.dffsOnTfItp.m2{t, 1}-min(rezG.dffsOnTfItp.m2{t, 1}),...
        rezR.dffsOnTfItp.m2{t, 1}-min(rezR.dffsOnTfItp.m2{t, 1}), 2000, 'normalized');
    
    %rs
    [rezG.dffsOnTfItp.rs{t, 1}, rezR.dffsOnTfItp.rs{t, 1}] = alignTwoDffTimeSeries(rezG.stimOnDffC.rs{t, 1}, rezR.stimOnDffC.rs{t, 1}, ...
        rezG.stimOnDffC.rs{t, 2}, rezR.stimOnDffC.rs{t, 2}, [-0.9 5]); 
    [rez.xcorrDffGR.rs{t, 1}, rez.xcorrDffGR.rs{t, 2}] = xcorr(rezG.dffsOnTfItp.rs{t, 1}-min(rezG.dffsOnTfItp.rs{t, 1}),...
        rezR.dffsOnTfItp.rs{t, 1}-min(rezR.dffsOnTfItp.rs{t, 1}), 2000, 'normalized');
    
    %ss
    [rezG.dffsOnTfItp.ss{t, 1}, rezR.dffsOnTfItp.ss{t, 1}] = alignTwoDffTimeSeries(rezG.stimOnDffC.ss{t, 1}, rezR.stimOnDffC.ss{t, 1}, ...
        rezG.stimOnDffC.ss{t, 2}, rezR.stimOnDffC.ss{t, 2}, [-0.9 5]); 
    [rez.xcorrDffGR.ss{t, 1}, rez.xcorrDffGR.ss{t, 2}] = xcorr(rezG.dffsOnTfItp.ss{t, 1}-min(rezG.dffsOnTfItp.ss{t, 1}),...
        rezR.dffsOnTfItp.ss{t, 1}-min(rezR.dffsOnTfItp.ss{t, 1}), 2000, 'normalized');
    
    %v1
    [rezG.dffsOnTfItp.v1{t, 1}, rezR.dffsOnTfItp.v1{t, 1}] = alignTwoDffTimeSeries(rezG.stimOnDffC.v1{t, 1}, rezR.stimOnDffC.v1{t, 1}, ...
        rezG.stimOnDffC.v1{t, 2}, rezR.stimOnDffC.v1{t, 2}, [-0.9 5]); 
    [rez.xcorrDffGR.v1{t, 1}, rez.xcorrDffGR.v1{t, 2}] = xcorr(rezG.dffsOnTfItp.v1{t, 1}-min(rezG.dffsOnTfItp.v1{t, 1}),...
        rezR.dffsOnTfItp.v1{t, 1}-min(rezR.dffsOnTfItp.v1{t, 1}), 2000, 'normalized');
    
    fprintf(sprintf("Completed trial#%d\n!", t))
end
