
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
function refTsI = check_timestamp_overlap(refTs, evtTs)
% Generate example data (replace with your actual data)
%refTs = faceCam_pt;
%evtTs = xStampSec_lick;

% Find the minimum and maximum timestamps of the first set
minRefTs = min(refTs);
maxRefTs = max(refTs);

refTsI = zeros(length(refTs), 1); 

% Select second timestamps that fall between the minimum and maximum of the first set
selectEvtTs = evtTs(evtTs(:, 1) >= minRefTs & evtTs(:, 1) <= maxRefTs, 1);

if ~isempty(selectEvtTs)
    for ts = 1:length(selectEvtTs)
        [~, minI] = min(abs(refTs - selectEvtTs(ts, 1)));
        refTsI(minI) = 1;
    end
end

end











