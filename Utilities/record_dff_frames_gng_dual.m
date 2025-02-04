function record_dff_frames_gng_dual(filePath, channel, varargin)
% filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_121324/task';
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.
% channel must be 'green' or 'red'

%% load tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

fileBeh = GrabFiles_sort_trials(strcat(channel, '_tbytDat_dff'), 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials(strcat(channel, '_tbytDat_dff'), 1, {filePath});
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
    if strcmpi(channel, 'green')
        pngSubDir = findFolderPattern(filePathTrials, ['dffG*', sprintf('trial_%d', t)]);
    elseif strcmpi(channel, 'red')
        pngSubDir = findFolderPattern(filePathTrials, ['dffR*', sprintf('trial_%d', t)]);
    end
    load(fullfile(pngSubDir, 'tbytDff.mat'), 'tbytDffsm');

    for i = 1:size(tbytDffsm,3)
        pngName = sprintf('frame_%d.png', i);
        figHandle = imageFrameWithNaNs(tbytDffsm(:, :, i), [-2 2]); hold on;

        % Turn off the axes
        axis off;
        % Set the figure's background to white
        set(gcf, 'Color', 'w');

        if isfield(tbytDat, 'frameStimI')
            if ~isempty(tbytDat(t).frameStimI) && tbytDat(t).frameStimI(i) % for auditory trials
                if tbytDat(t).evtType == 1 % common auditory stim (draw 45 degree lines at the upper left corner)
                    text(2, 3, 'Tone 1', 'Color', 'white', 'FontSize', 16);
                elseif tbytDat(t).evtType == 2 % common auditory stim (draw 135 degree lines at the upper left corner)
                    text(2, 3, 'Tone 2', 'Color', 'white', 'FontSize', 16);
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fullPath = findFolderPattern(baseFolder, pattern)
% FIND FOLDER in a given directory that matches a given pattern exactly.
%
% Inputs:
%   baseFolder - The parent directory to search within
%   pattern    - A string with '*' wildcards (e.g., 'dffG*trial_1') or
%                a cell array of substrings (e.g., {'dffG', 'trial_1'}).
%
% Output:
%   fullPath - Full path to the matched folder or empty if not found.

% List all subdirectories
dirs = dir(baseFolder);

% Filter out non-directories
dirs = dirs([dirs.isdir]);

% Remove '.' and '..' from results
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));

% Extract folder names
folderNames = {dirs.name};

if ischar(pattern)
    % Convert wildcard '*' to regex '.*'
    regexPattern = strrep(pattern, '*', '.*');

    % Ensure exact number matching for digits (avoid partial matches like 'trial_11')
    regexPattern = regexprep(regexPattern, '(\d+)', '$1(?!\d)');

    % Perform regex match
    matchIdx = ~cellfun(@isempty, regexp(folderNames, ['^' regexPattern '$'], 'once'));

elseif iscell(pattern)
    % Match substrings with exact number matching
    matchIdx = cellfun(@(name) all(contains(name, pattern)) && ...
        all(cellfun(@(p) ~contains(name, [p, '0':'9']), pattern)), folderNames);

else
    error('Pattern must be a string (wildcard) or a cell array of substrings.');
end

% Return the matched folder path
if any(matchIdx)
    fullPath = fullfile(baseFolder, folderNames{matchIdx});
else
    fullPath = ''; % Return empty if no match found
end
end



