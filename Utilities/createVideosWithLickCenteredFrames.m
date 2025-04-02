%This function creates videos comprising frames centered around licking and
% other orofacial behavioral events of interest  

filePaths = {'/Volumes/buschman/Rodent Data/dualImaging_parkj/m1237_GCAMP/m1237_100124/task', ...           % cropped 
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1237_GCAMP/m1237_101024/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO/m1092_100924/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO/m1092_101024/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1094_jRGECO/m1094_100924/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1044_jRGECO_GRABda/m1044_121924/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1044_jRGECO_GRABda/m1044_122024/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122024/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122324/task', ...             
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1048_jRGECO_GRABda/m1048_122324/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_122024/task', ...
             }; 

fileSavePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectVideo/orofacialToTrainDLC'; 

% User-defined parameters
win = [-1 2]; % time window for the new edited videos  

for f = 1:length(filePaths)
    % Grab behavioral data
    tbytDatPath = GrabFiles_sort_trials('_tbytDat_parseGng', 0, fullfile(filePaths(f), 'Matfiles')); 
    load(tbytDatPath{1}, 'tbytDat')

    filePath_evt = findFileInSubfolders(filePaths{f}, 'evtInS.mat');
    load(fullfile(filePath_evt), 'evtInS');
    
    header = extract_date_animalID_header(filePaths{f}); 


    % Specify trials of interest
    hitLickN = cellfun(@numel, {tbytDat.hitLicks});  
    [~, hitLickI] = sort(hitLickN); 
    
    trialId = hitLickI(end); 
        
    if ~isempty(tbytDat(trialId).hitLicks)
       % Identify faceCam and frames 
       tBound = tbytDat(trialId).hitLicks(1)+win; 
       faceCamTrialTs = tbytDat(trialId).faceCam-tbytDat(trialId).evtOn;
       faceCamTrialTsI = faceCamTrialTs >= tBound(1) & faceCamTrialTs <= tBound(end); 

       faceCamI = tbytDat(trialId).topRedExpTrainI;
       faceCamTs = evtInS.faceCam(evtInS.faceCam(:, 2)==faceCamI, 1);
       correspondingFrameI = interp1(faceCamTs, 1:length(faceCamTs), tbytDat(trialId).faceCam, 'nearest'); % tbytDat(t).resampledVidFrameTs
    
       faceCamFrameId = correspondingFrameI(faceCamTrialTsI); 

       filePath_ofVid = GrabFiles_sort_trials('_cropped_orofacial', 0, filePaths(f)); 
       fileList_ofVid = GrabFiles_sort_trials('_cropped_orofacial.mp4', 0, filePath_ofVid); 
       fileList_ofVid{faceCamI}; 

       fileSavePathFull = fullfile(fileSavePath, [header sprintf('_trial%d_hitLickAligned', trialId)]); 

       writeVideoSelectFrames(fileList_ofVid{faceCamI}, faceCamFrameId, fileSavePathFull)

    end

end