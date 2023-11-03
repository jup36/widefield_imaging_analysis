function record_dff_frames_gng(filePath, varargin)
% filePath = 'Z:\Rodent Data\Behavioral_dynamics_cj\DA008\DA008_101823';
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.

%% load tbytDat
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

%% sort trials to work with
if isempty(varargin)
    trials = 1:length(tbytDat);
else
    trials = [varargin{1}(:)]';
end

%%
filePathTrials = fullfile(filePath, strcat(header, '_trials'));

%% record frames first
close all;
for t = trials % 1:length(tbytDat) % trial
    pngSubDir = fullfile(filePathTrials, sprintf('trial_%d', t));
    load(fullfile(pngSubDir, 'tbytDff.mat'), 'tbytDffsm');

    for i = 1:size(tbytDffsm,3)
        pngName = sprintf('frame_%d.png', i);
        figHandle = imageFrameWithNaNs(tbytDffsm(:, :, i), [-2 2]); hold on;

        % Turn off the axes
        axis off;
        % Set the figure's background to white
        set(gcf, 'Color', 'w');

        if isfield(tbytDat, 'frameStimI')
            if ~isempty(tbytDat(t).frameStimI) && tbytDat(t).frameStimI(i) % for visual trials
                if tbytDat(t).rewardTrI == 1 % common visual stim (draw 0 degree lines at the upper left corner)
                    insertgrating0(figHandle, tbytDffsm(:, :, i))
                elseif tbytDat(t).punishTrI == 1 % common visual stim (draw 90 degree lines at the upper left corner)
                    insertgrating90(figHandle, tbytDffsm(:, :, i))
                end
                hold off;
            end
        end
        % Capture the current figure with dots
        frameLabled = getframe(gca);
        frameLabled = frameLabled.cdata;

        % save the figure with dots

        imwrite(frameLabled, fullfile(pngSubDir, pngName));

        fprintf('Frame #%d of trial #%d is labeled and saved.\n', i, t);
        close all
    end
end
end


