%% ==========================
% Setup: load rezSubspaceSim
% ===========================
load(compatiblepath("Z:\Rodent Data\dualImaging_parkj\collectData\glmTdrSubspaceRezCollection_031826.mat"), ...
    "rezSubspaceSim")

%% =========================
% Setup: mouse IDs, learning groups, and plotting colors
% =========================

% Extract mouse IDs from the per-mouse results structure
mIdC = cellfun(@(a) a.mouseId, rezSubspaceSim.perMouse, 'UniformOutput', false);

% Define fast- and slow-learning animal groups
fastC = {'m1044', 'm1045', 'm1092', 'm1094'};
slowC = {'m1048', 'm1049', 'm1613', 'm1859', 'm1873'};

% Convert group membership into logical indices matching rezSubspaceSim.perMouse
fastIDs = ismember(mIdC, fastC);
slowIDs = ismember(mIdC, slowC);

% Define group colors: row 1 = fast, row 2 = slow
colorMatFS = [157, 0, 255; 80, 200, 120] ./ 255;

%% =========================
% Plot and save overlap-to-reference trajectories for multiple subspaces
% =========================

% List of subspaces to visualize
subspaceNames = {'GoToneOn', 'NoGoToneOn', 'ToneOffGo', 'ToneOffNoGo', 'GoTone', 'NoGoTone'};

% Optional: store output handles from each plot
hSubspace = struct();

for iS = 1:numel(subspaceNames)

    % Current subspace name
    subspaceName = subspaceNames{iS};

    % Extract overlap trajectories for the current subspace
    overlapC = cellfun( ...
        @(a) a.relativeToExpert.(subspaceName).overlap, ...
        rezSubspaceSim.perMouse, ...
        'UniformOutput', false);

    % Plot and save
    hSubspace.(subspaceName) = subspaceOverlapRelToRef_fastSlow( ...
        overlapC, fastIDs, slowIDs, colorMatFS, ...
        'ylim', [0.2 0.9], ...
        'axisTight', true,...
        'MarkerSizeGroup', 100, ...
        'TitleStr', sprintf('%s overlap', subspaceName), ...
        'printFigLogic', true, ...
        'figSaveName', sprintf('refSession_%s_overlap_fastSlow', subspaceName));
end

%% =========================
% Cross-trial-type pair names
% =========================
pairNames = rezSubspaceSim.perMouse{find(~cellfun(@isempty, rezSubspaceSim.perMouse),1)} ...
    .withinSession.crossTrialType.pairNames;

%% =========================
% Plot each pair across sessions
% =========================
hWithin = struct();

for iP = 1:numel(pairNames)

    pairName = pairNames{iP};

    % Extract one overlap vector per mouse for the current pair
    pairOverlapC = cellfun(@(a) a.withinSession.crossTrialType.overlap(:,iP), ...
        rezSubspaceSim.perMouse, 'UniformOutput', false);

    % Plot and save
    hWithin.(matlab.lang.makeValidName(pairName)) = withinSessionOverlap_fastSlow( ...
        pairOverlapC, fastIDs, slowIDs, colorMatFS, ...
        'leftAlignFirstValid', true, ...
        'ylim', [0.2 0.95], ...
        'axisTight', true, ...
        'MarkerSizeGroup', 100, ...
        'TitleStr', strrep(pairName, '__vs__', ' vs '), ...
        'printFigLogic', true, ...
        'figSaveName', sprintf('withinSession_%s_fastSlow', pairName));
end