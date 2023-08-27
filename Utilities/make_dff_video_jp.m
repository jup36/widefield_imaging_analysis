
%% whereabouts
filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623';
filePath_nidq = GrabFiles_sort_trials('_g', 0, {filePath}); 
if isempty(filePath_nidq)
    filePath_nidq = uigetdir(filePath, 'Select a folder');
end
[file_list_first_stack,folder_list_raw] = GrabFiles_sort_trials('ome_stack.mat',1, ... % use GrabFiles_sort_trials to sort both files and folders 
    {'/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623/DA001_072623_img'}, true);
trNumb = 13; 
filePath_png = fullfile(strcat(folder_list_raw{trNumb},'/png'));
if exist(filePath_png, 'dir') == 0
    mkdir(filePath_png);
end

filePath_vid = fullfile(strcat(folder_list_raw{trNumb},'/vid'));
if exist(filePath_vid, 'dir') == 0
    mkdir(filePath_vid);
end

%% align frames to task timestamps
evtInS = timestamp_behav_events(filePath_nidq{1}, true, 'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'photoDiode'); % behavioral events

% load trial-by-trial behav data
%[~, mainFolder] = fileparts(filePath); 
%[fileBeh, ~] = GrabFiles_sort_trials(mainFolder(2:end),0, ... % use GrabFiles_sort_trials to sort both files and folders 
%    {'/Volumes/buschman/Users/Caroline/NADA_dynamics/data'}, true);
%tbytDat = load(fullfile(fileBeh{1}), 'data');
%tbytDat = tbytDat.('data');

% load a stack
load(file_list_first_stack{trNumb}, 'stack'); 

% determine the frame order
stack_odd = stack(:, :, 1:2:end); 
stack_even = stack(:, :, 2:2:end); 

frameTimeN_total = evtInS.cmosExp(evtInS.cmosExp(:, 2)==trNumb); 

if mean(stack_odd(:)) >= mean(stack_even(:))
    stack_b = stack_odd; 
    frameTimeN = frameTimeN_total(1:2:end);  
else
    stack_b = stack_even;
    frameTimeN = frameTimeN_total(2:2:end);  
end

% map temporal events to camera pulses
frameLickI = check_timestamp_overlap(frameTimeN, evtInS.lick);
frameWaterI = check_timestamp_overlap(frameTimeN, evtInS.water);
%frameAirpuffI = check_timestamp_overlap(frameTimeN, evtInS.airpuff);

% timestamp each frame relative to the stim onset
frameStimOnI = check_timestamp_overlap(frameTimeN, evtInS.photoDiode(trNumb, 1));
frameStimOffI = check_timestamp_overlap(frameTimeN, evtInS.photoDiode(trNumb, 2));
frameStimI = zeros(length(frameStimOnI), 1);
frameStimI(find(frameStimOnI, 1):find(frameStimOffI, 1), 1) = 1;
reframeTimeN = round(frameTimeN-evtInS.photoDiode(trNumb, 1), 2);

%basic dff
baseI = check_timestamp_overlap(frameTimeN, evtInS.photoDiode(trNumb, 1)-2);
avg_b = nanmean(stack_b(:,:,1:find(baseI)),3); % baseline average (0.5 sec)
dff = (double(stack_b)-avg_b)./avg_b; 

%% examine the data
%temp = reshape(dff,size(dff,1)*size(dff,2),size(dff,3));
%imagesc(temp)

%histrogram
%figure; hold on; 
%histogram(dff(:))

%% record frames first 
textPos = {[5, 63], [33, 63], [60, 63]};
close all;
for i = 1:size(dff,3)
    imagesc(dff(:,:,i),[-0.04, 0.04]); hold on; 
    %set(gca, 'XTickLabel', '', 'YTickLabel', '');
    %title(sprintf('frame %d of %d',i,size(dff,3)));
    
    if frameStimI(i)
        % draw the grating pattern
        rectangle('Position', [[5, 5], [7 1]], 'EdgeColor', 'none', 'FaceColor', 'white');
        rectangle('Position', [[5, 6], [7 1]], 'EdgeColor', 'none', 'FaceColor', 'black');
        rectangle('Position', [[5, 7], [7 1]], 'EdgeColor', 'none', 'FaceColor', 'white');
        rectangle('Position', [[5, 8], [7 1]], 'EdgeColor', 'none', 'FaceColor', 'black');
        rectangle('Position', [[5, 9], [7 1]], 'EdgeColor', 'none', 'FaceColor', 'white');
    end
    hold off; 
    
    % Capture the current figure with dots
    frameLabled = getframe(gcf);
    frameLabled = frameLabled.cdata;

    % save the figure with dots
    png_name = sprintf('frame_%d.png', i);
    imwrite(frameLabled, fullfile(filePath_png, png_name));

    fprintf('Frame #%d is labeled and saved.\n', i);
end

%% write video
frameRate = 10; 
labeled_videos = findAsManyFilesWithString(fullfile(filePath_vid), '_labeled'); 
vidname.tr = sprintf('_tr%d', trNumb); 
vidname.order = sprintf('_%d', length(labeled_videos)+1);
vidname.fr = sprintf('_FR%d', frameRate); 
vidname.full = fullfile(filePath_vid, strcat(mainFolder, vidname.tr, vidname.order, vidname.fr, '.mp4')); 

labeledVideo = VideoWriter(vidname.full, 'MPEG-4');
labeledVideo.FrameRate = frameRate;
playspeed = labeledVideo.FrameRate/(size(dff,3)/15); 
open(labeledVideo);

pngFiles = GrabFiles_sort_trials('frame_', 0, {filePath_png}); 

for ff = 1:length(pngFiles)
    % Load PNG image;
    img = imread(pngFiles{ff}); 

    % Add playspeed to the image
    playSpeedStr = sprintf('playspeed x%.1f', playspeed);
    imgWithText = insertText(img, [165, 640], playSpeedStr, 'FontSize', 24, 'TextColor', 'white', 'BoxColor', [1, 0, 1], 'BoxOpacity', 0.8);

    % Add frametime to the image
    frameTimeStr = sprintf('%.2fs', reframeTimeN(ff));
    imgWithText = insertText(imgWithText, [165, 690], frameTimeStr, 'FontSize', 24, 'TextColor', 'white');

    if frameLickI(ff)
        imgWithText = insertText(imgWithText, [265, 690], 'Lick', 'FontSize', 24, 'TextColor', 'white');
    end

    if frameWaterI(ff)
        imgWithText = insertText(imgWithText, [365, 690], 'Water', 'FontSize', 24, 'TextColor', 'white');
    end
    
    frame = im2frame(imgWithText); 

    % write the labeled frame to the output video
    writeVideo(labeledVideo, frame);
    close; 
    fprintf('Frame #%d is loaded and written.\n', ff);

end

% close the output video
close(labeledVideo);

disp('Video creation complete.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evtInS = timestamp_behav_events(path_raw, redetect_logic, varargin)
%This function takes the number of behavioral variables as varargin and
% gets the timestamps (in sec) of the variables and saves them as a
% structure named evtInS in the files origianl nidq directory.
% TO DO: Add more variables in the SWITCH CASE block.

% Load bin file
binFile = dir(fullfile(path_raw, '*.nidq.bin')); % look for nidq.bin file
binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
binFile = binFile(binFileI);

if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Load evtInS
if redetect_logic
    evtIns = struct;
else
    if ~isempty(findFileWithString(path_raw, 'evtInS'))
        load(fullfile(path_raw, 'evtInS'))
        evtInS_fields = fieldnames(evtInS);
        evtInS_isfield = cell2mat(cellfun(@(a) any(strcmp(a, evtInS_fields)), varargin, 'un', 0));
        if sum(evtInS_isfield)==length(varargin)
            return; % if all the events were already detected and available just return
        else % in case there's still to be detected and appended
            varargin = varargin(~evtInS_isfield);
        end
    else
        evtInS = struct;
    end
end

% Parse the corresponding metafile
meta  = ReadMeta(binName, path_raw); % get the meta data (structure)
sample_rate=str2double(meta.niSampRate);
nSamp = floor(SampRate(meta));          % sampling rate (default: 25kHz)
%totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds
%channels = textscan(meta.acqMnMaXaDw,'%n %n %n %n','Delimiter',',');

output = struct;

%% Main (load raw traces and perform event detection)
if exist(fullfile(path_raw,'gainCorrectRawTraces.mat'), 'file')~=2
    p = parse_input_PP(nidq_data_file, {});
    behaviorTimestampsPP(p);
end

if exist(fullfile(path_raw,'gainCorrectRawTraces.mat'), 'file')==2
    for jj = 1:length(varargin)
        fieldName = sprintf('variable_%d', jj);
        output.(fieldName) = load(fullfile(path_raw,'gainCorrectRawTraces.mat'), varargin{jj});
        assert(isfield(output.(fieldName), varargin{jj})) % ensure the raw trace is loaded properly

        switch varargin{jj}
            case 'lick'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'lick'))
                    lick = output.(fieldName).lick;
                    timeStamp = (1:length(lick))';
                    timeStampSec = timeStamp./sample_rate;
                    lickIdx = detecteventbythreshold(lick, nSamp, 50, 'stdFactor', 3, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true);
                    evtInS.lick = timeStampSec(lickIdx);
                    fprintf('completed lick detection!\n');
                end
            case 'faceCam'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'faceCam'))
                    faceCam = output.(fieldName).faceCam;
                    timeStamp = (1:length(faceCam))';
                    timeStampSec = timeStamp./sample_rate;
                    [faceCamRiseIdx, ~, faceCamTrainIdx] = detecteventbythreshold_noAbs(faceCam, nSamp, 4, 'stdFactor', 2.5, 'plotRez', false, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', false); % camera trigger
                    evtInS.faceCam = [timeStampSec(faceCamRiseIdx), faceCamTrainIdx'];
                    fprintf('completed faceCam detection!\n');
                end

            case 'cmosExp'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'cmosExp'))
                    cmosExp = output.(fieldName).cmosExp;
                    timeStamp = (1:length(cmosExp))';
                    timeStampSec = timeStamp./sample_rate;
                    [cmosExpRiseIdx, ~, cmosExpTrainIdx] = detecteventbythreshold_noAbs(cmosExp, nSamp, 20, 'stdFactor', 1.5, 'plotRez', false, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', false); % CMOS trigger
                    evtInS.cmosExp = [timeStampSec(cmosExpRiseIdx), cmosExpTrainIdx'];
                    fprintf('completed cmosExp detection!\n');
                end
            case 'water'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'water'))
                    water = output.(fieldName).water;
                    timeStamp = (1:length(water))';
                    timeStampSec = timeStamp./sample_rate;
                    waterIdx = detecteventbythreshold_noAbs(water, nSamp, 50, 'stdFactor', 1.5, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % CMOS trigger
                    evtInS.water = timeStampSec(waterIdx);
                    fprintf('completed water detection!\n');
                end

            case 'photoDiode'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'photoDiode'))
                    photoDiode = output.(fieldName).photoDiode;
                    timeStamp = (1:length(photoDiode))';
                    timeStampSec = timeStamp./sample_rate;
                    [photoDiodeRiseIdx, photoDiodeFallIdx ] = detecteventbythreshold_noAbs(photoDiode, nSamp, 50, 'stdFactor', 1.5, 'plotRez', true, ...
                        'chunkPulses', false, 'correctLongPulse', true, 'cutoffShort', true, 'short', 1, 'findEarlyOnset', true, 'earlyOnsetWindow', 0.5); % CMOS trigger
                    evtInS.photoDiode = [timeStampSec(photoDiodeRiseIdx), timeStampSec(photoDiodeFallIdx)];
                    fprintf('completed photoDiode detection!\n');
                end

            case 'airpuff'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'airpuff'))
                    airpuff = output.(fieldName).airpuff;
                    if isempty(airpuff)
                        evtInS.airpuff = []; 
                    else
                        timeStamp = (1:length(airpuff))';
                        timeStampSec = timeStamp./sample_rate;
                        airpuffIdx = detecteventbythreshold_noAbs(airpuff, nSamp, 50, 'stdFactor', 1.5, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % CMOS trigger
                        evtInS.airpuff = timeStampSec(airpuffIdx);
                        fprintf('completed airpuff detection!\n');
                    end
                end
        end
    end
end

%% save evtInS
if exist(fullfile(path_raw, "evtInS"), 'file') == 2
    save(fullfile(path_raw, "evtInS"), 'evtInS', '-append')
else
    save(fullfile(path_raw, "evtInS"), "evtInS")
end

fprintf('completed event detection and saved events in the nidq folder!\n');

end









