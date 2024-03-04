function make_dff_videos_auditory(tbytDat, filePathTrials, frameRate, varargin)
close all;
filePathTrials = GrabFiles_sort_trials('trial', 0, {filePathTrials}); % use GrabFiles_sort_trials to sort both files and folders

if isempty(varargin)
    trials = 1:length(filePathTrials);
else
    % To match the trial number use this!
    tN = cellfun(@(x) str2double(regexp(x, 'trial_(\d{1,3})', 'tokens', 'once')), filePathTrials);
    trials = [varargin{1}(:)]';
    filePathTrials = filePathTrials(ismember(tN, trials));
end

if isempty(filePathTrials)
    warning('There is no trials ready to be filmed!')
end

% To match the trial number use this!
%tN = cellfun(@(x) str2double(regexp(x, 'trial_(\d{1,3})', 'tokens', 'once')), filePathTrials);
%t = find(tN == 622);

for t = 1:length(filePathTrials) % trials

    [~, vidName] = fileparts(filePathTrials{t});

    % Extract the number after 'trial_'
    matchedNum = regexp(vidName, 'trial_(\d+)', 'tokens');
    % Convert the matched number from cell to double
    trialNum = str2double(matchedNum{1}{1});
    frameTrel = tbytDat(trialNum).frameTrel; % get frame times relative to evt

    writeDffVid_simple(filePathTrials{t}, vidName, frameRate, frameTrel)

    fprintf('Video #%d of total %d videos is written.\n', t, length(filePathTrials));
end




    function writeDffVid_timeWin(pngPath, vidName, frameRate, frameTrel, timeWin)
        pngs = GrabFiles_sort_trials('frame_', 0, {pngPath});
        labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
        labeledVid.FrameRate = frameRate;
        playspeed = labeledVid.FrameRate/15;
        open(labeledVid);

        frameTrelWin = frameTrel >= timeWin(1) & frameTrel <= timeWin(2); 

        for f = 1:length(pngs) % frames
            % load png image
            if frameTrelWin(f)
            img = imread(pngs{f});

            % Add playspeed to the image
            playSpeedStr = sprintf('playspeed x%.1f', playspeed);
            imgWithText = insertText(img, [10, 607], playSpeedStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [1, 0, 1], 'BoxOpacity', 0.8);

            % Add frametime to the image
            frameTimeStr = sprintf('%.2fs', frameTrel(f));
            imgWithText = insertText(imgWithText, [700, 607], frameTimeStr, 'FontSize', 22, 'TextColor', 'white');

            imshow(imgWithText);

            frame = im2frame(imgWithText);

            % write the labeled frame to the output video
            writeVideo(labeledVid, frame);
            %fprintf('Frame #%d is loaded and written.\n', f);
        
            end
        end
        close(labeledVid);
    end

    function writeDffVid_simple(pngPath, vidName, frameRate, frameTrel)
        pngs = GrabFiles_sort_trials('frame_', 0, {pngPath});
        labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
        labeledVid.FrameRate = frameRate;
        playspeed = labeledVid.FrameRate/15;
        open(labeledVid);

        for f = 1:length(pngs) % frames
            % load png image
            img = imread(pngs{f});

            % Add playspeed to the image
            playSpeedStr = sprintf('playspeed x%.1f', playspeed);
            imgWithText = insertText(img, [10, 607], playSpeedStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [1, 0, 1], 'BoxOpacity', 0.8);

            % Add frametime to the image
            frameTimeStr = sprintf('%.2fs', frameTrel(f));
            imgWithText = insertText(imgWithText, [700, 607], frameTimeStr, 'FontSize', 22, 'TextColor', 'white');

            imshow(imgWithText);

            frame = im2frame(imgWithText);

            % write the labeled frame to the output video
            writeVideo(labeledVid, frame);
            %fprintf('Frame #%d is loaded and written.\n', f);
        end
        close(labeledVid);
    end
end
