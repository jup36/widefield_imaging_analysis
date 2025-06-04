function frameTimes_for_dffFaceVideos(filePath)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224';
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.
% Note that face videos

%% load tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

% filePathTrials = GrabFiles_sort_trials('_trials', 0, {filePath});

fileBehG= GrabFiles_sort_trials('green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehG{1})
    fileBehG = GrabFiles_sort_trials('green_tbytDat_dff', 1, {filePath});
end
tbytDatG = load(fullfile(fileBehG{1}), 'tbytDat');
tbytDatG = tbytDatG.('tbytDat'); clearvars tbytDat

fileBehR= GrabFiles_sort_trials('red_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehR{1})
    fileBehR = GrabFiles_sort_trials('red_tbytDat_dff', 1, {filePath});
end
tbytDatR = load(fullfile(fileBehR{1}), 'tbytDat');
tbytDatR = tbytDatR.('tbytDat'); clearvars tbytDat

tbytDatG = addRewardTrI(tbytDatG);
tbytDatR = addRewardTrI(tbytDatR);
tbytDatG = addPunishTrI(tbytDatG);
tbytDatR = addPunishTrI(tbytDatR);

% load evtInS
fileNidq = GrabFiles_sort_trials('_g', 0, {filePath});
if ~isempty(fileNidq) && isscalar(fileNidq)
    fileEvt = GrabFiles_sort_trials('evtInS', 0, fileNidq(1));
    load(fullfile(fileEvt{1}), 'evtInS');
else
    [evtInS_file, evtInS_folder] = uigetfile('*.mat', 'Select the evtInS file for faceCam!', filePath);
    load(fullfile(evtInS_folder, evtInS_file), 'evtInS')
end

%% list face videos
faceVidPath = GrabFiles_sort_trials([header, '*', '_vid'], 0, {filePath});
faceVids = GrabFiles_sort_trials('behvid', 0, faceVidPath);

%% record frames first
downSampleFactor = 4; % default downsample factor to downsample the faceCam frames
close all;
for t = 1:length(tbytDatG)
    % interpolate dff (10Hz) to match the faceCam sampling rate (200Hz)
    assert(tbytDatG(t).sideGreenExpTrainI==tbytDatR(t).topRedExpTrainI)
    faceCamI = tbytDatG(t).sideGreenExpTrainI;
    faceCamTs = evtInS.faceCam(evtInS.faceCam(:, 2)==faceCamI, :);

    targetFaceTsFstI = find(faceCamTs(:, 1) < min(tbytDatG(t).frameT(1), tbytDatR(t).frameT(1)), 1, 'last');
    targetFaceTsLastI = find(faceCamTs(:, 1) > max(tbytDatG(t).frameT(end), tbytDatR(t).frameT(end)), 1, 'first');
    faceFramesToRead = targetFaceTsFstI:downSampleFactor:targetFaceTsLastI; % downsampled

    v = VideoReader(faceVids{faceCamI});
    totalFaceFrames = v.Duration*v.FrameRate;

    if ~isempty(targetFaceTsFstI) && ~isempty(targetFaceTsLastI) && targetFaceTsLastI <= totalFaceFrames
        targetFaceTs = faceCamTs(faceFramesToRead, 1); % downsampled
        tbytDatG(t).resampledVidFrameTs = targetFaceTs;
        tbytDatR(t).resampledVidFrameTs = targetFaceTs;
    end
    clear v;
    fprintf('Trial #%d is completed.\n', t);
end

save(fullfile(fileBehG{1}), 'tbytDatG', '-append') % to save tbytDat(t).faceCamTrelFrames
save(fullfile(fileBehR{1}), 'tbytDatR', '-append') % to save tbytDat(t).faceCamTrelFrames

end