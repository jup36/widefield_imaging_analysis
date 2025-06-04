function dffPostprocess_auditory_gng_averaging_dual_6OHDA(filePath, channel)

%% Load Data
header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match', 'once');

% Load behavior data based on channel
channel_suffix = ['_' channel '_tbytDat_dff'];
fileBeh = GrabFiles_sort_trials(channel_suffix, 0, {fullfile(filePath, 'Matfiles')});
load(fullfile(fileBeh{1}), 'tbytDat');

% Load event data
evt_file = GrabFiles_sort_trials('evtInS', 0, GrabFiles_sort_trials('_g', 0, {filePath}));
if isempty(evt_file)
    [evt_file, evt_folder] = uigetfile('*.mat', 'Select evtInS file', filePath);
    evt_file = fullfile(evt_folder, evt_file);
else
    evt_file = fullfile(evt_file{1});
end
load(evt_file, 'evtInS');

% Organize manual water timestamps
tbytDat = organizeManualWaterIntoTbytDat(evtInS.manualWater, tbytDat);
manualWaterC = {tbytDat.manualW};

mWaterV = NaN(1, numel(tbytDat));
hasMwater = ~cellfun('isempty', {tbytDat.manualW});
mWaterV(hasMwater) = cellfun(@(a) a(1), {tbytDat(hasMwater).manualW});

hitLickV = NaN(1, numel(tbytDat));
hasHit = ~cellfun('isempty', {tbytDat.hitLicks});
hitLickV(hasHit) = cellfun(@(a) a(1), {tbytDat(hasHit).hitLicks});

voluntaryHitLickI = isnan(mWaterV) & ~isnan(hitLickV) | hitLickV<mWaterV; 

% Pre-allocate with NaN
hitLickC = NaN(1, numel(tbytDat));

% Logical mask for the trials that actually have hit-licks
hasHit = ~cellfun('isempty', {tbytDat.hitLicks});

% Fill only those indices
hitLickC(hasHit) = cellfun(@(a) a(1), {tbytDat(hasHit).hitLicks});

% Load preprocessed dffs
load(fullfile(filePath, 'Matfiles', sprintf('%s_%s_dff_smCollect.mat', header, channel)), 'dffsmCell');

% Load Allen maps and masks
load('allenDorsalMap.mat', 'dorsalMaps', 'motorMask', 'smotorMask', 'ssMask', 'vMask', 'rsMask');
masks = {motorMask, smotorMask, ssMask, vMask, rsMask};
regions = {'m1', 'm2', 'ss', 'v1', 'rs'};

% Generate hemispheric masks correctly
hemiMasksL = cellfun(@(m) m & [true(size(m,1), floor(size(m,2)/2)), false(size(m,1), ceil(size(m,2)/2))], masks, 'UniformOutput', false);
hemiMasksR = cellfun(@(m) m & [false(size(m,1), floor(size(m,2)/2)), true(size(m,1), ceil(size(m,2)/2))], masks, 'UniformOutput', false);

% Load Allen transformation object
img_folder = GrabFiles_sort_trials('_img', 0, {filePath});
transParams_file = find_keyword_containing_files(img_folder{1}, ['transParamsAllen_' channel], 'recursive', true);
load(fullfile(transParams_file{1}), 'transParams');

%% Main Processing Loop
% Pre-Initialize rez structure
nTrials      = numel(tbytDat);
eventFields  = {'stimOnDffC','postStimLickDff','firstWater', ...
                'hitLickFst','hitLickLast','hitLickToneAlign','consumeLick'};
regionLabels = strcat(regions, 'L');           %  m1L … rsL
regionLabels = [regionLabels, strcat(regions,'R')];  % add m1R … rsR

rez = struct;
for ef = 1:numel(eventFields)
    for rl = 1:numel(regionLabels)
        rez.(eventFields{ef}).(regionLabels{rl}) = cell(nTrials,2);   % pre-fill with []
    end
end

for t = 1:length(tbytDat)

    % Align to AllenCCF
    dff_aligned = alignStackToAllenKabsch(dffsmCell{t}, dorsalMaps.dorsalMap, transParams.tformObj);

    % Apply region masks
    for k = 1:numel(regions)
        dffL{k} = apply2DMaskTo3DStack(dff_aligned, hemiMasksL{k});
        dffR{k} = apply2DMaskTo3DStack(dff_aligned, hemiMasksR{k});
    end

    frT = tbytDat(t).frameTrel;

    %% Alignment blocks
    % Stim alignment
    rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'stimOnDffC', 0, frT, [-0.9 5], t);

    % Post-stim licks
    for jj = 1:length(tbytDat(t).postStimChunk)
        evt_time = tbytDat(t).postStimChunk{jj}(1);
        if evt_time + 1 < frT(end)
            rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'postStimLickDff', evt_time, frT, [-1 1], t, jj);
        end
    end

    % Align to water
    if ~isempty(tbytDat(t).water)
        firstWater = tbytDat(t).water(1) - tbytDat(t).evtOn;
        rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'firstWater', firstWater, frT, [-1 1], t);
    end

    % Align to hit licks
    if ~isempty(tbytDat(t).hitLicks)
        rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'hitLickFst', tbytDat(t).hitLicks(1), frT, [-1 1], t);
        rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'hitLickLast', tbytDat(t).hitLicks(end), frT, [-1 1], t);
        rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'hitLickToneAlign', 0, frT, [-0.9 1], t);

        if ~isempty(tbytDat(t).consumeLicks)
            rez = alignAndStoreEvent(dffL, dffR, regions, rez, 'consumeLick', tbytDat(t).consumeLicks(1), frT, [-1 1], t);
        end
    end
    fprintf("processed trial #%d\n", t);
end

%% Save Results
save(fullfile(filePath, 'Matfiles', sprintf('%s_%s_dff_evtAligned_regionMask.mat', header, channel)), 'rez', 'voluntaryHitLickI');

end

%% helper function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rez = alignAndStoreEvent(dffL, dffR, regions, rez, eventField, eventTime, frT, window, trialIdx, varargin)

if ~isfield(rez, eventField)
    rez.(eventField) = struct();
end

for k = 1:numel(regions)
    regionL = [regions{k}, 'L'];
    regionR = [regions{k}, 'R'];

    [alignedL, alignedTL] = alignToEvent(dffL{k}, eventTime, frT, window);
    [alignedR, alignedTR] = alignToEvent(dffR{k}, eventTime, frT, window);

    if ~isempty(varargin)
        jj = varargin{1};
        rez.(eventField).(regionL){trialIdx, 1}{jj, 1} = alignedL;
        rez.(eventField).(regionL){trialIdx, 1}{jj, 2} = alignedTL;
        rez.(eventField).(regionR){trialIdx, 1}{jj, 1} = alignedR;
        rez.(eventField).(regionR){trialIdx, 1}{jj, 2} = alignedTR;
    else
        rez.(eventField).(regionL){trialIdx, 1} = alignedL;
        rez.(eventField).(regionL){trialIdx, 2} = alignedTL;
        rez.(eventField).(regionR){trialIdx, 1} = alignedR;
        rez.(eventField).(regionR){trialIdx, 2} = alignedTR;
    end
end
end

function [evtAlignedDffTs, evtAlignedTs] = alignToEvent(dffTs, eventTimeToAlign, frameT, timeWin)
% NOTE: Modified the function to return outcomes even if the timeWin goes
% out of bounds of the frame time! (2/4/2025, Junchol Park)

win = timeWin + eventTimeToAlign;

frameI = frameT >= min(win) & frameT <= max(win);
evtAlignedDffTs = dffTs(frameI);
evtAlignedDffTs = evtAlignedDffTs(:).'; % make it a row vector
evtAlignedTs = frameT(frameI)-eventTimeToAlign;
evtAlignedTs = evtAlignedTs(:).'; % make it a row vector

end


function [imagStackMaskedMean, imgStackMasked] = apply2DMaskTo3DStack(imgStack, mask2D)
% Replicate the 2D mask to match the 3D stack dimensions
replicatedMask = repmat(mask2D, [1, 1, size(imgStack, 3)]);

% Apply the mask: replace values in the stack where the mask is zero (or false) with NaN
imgStack(~replicatedMask) = NaN;

imgStackMasked = imgStack;

imagStackMaskedMean = squeeze(nanmean(nanmean(imgStackMasked, 1), 2));

end




