function dffPostprocess_auditory_gng_averaging_dual_regionMask(filePath, channel)
%--------------------------------------------------------------------------
%  Align dF/F stacks to Allen CCF, extract regional traces, and align them
%  to behavioural events for an auditory go/no-go task.
%
%  *** 2025-06-12 update ***
%  – Removed left/right hemisphere split.  Each region is treated as a
%    single mask (m1, m2, ss, v1, rs).
%--------------------------------------------------------------------------

%% 1. Load behaviour and event data
header = regexp(filePath,'m\d{1,4}_\d{6}','match','once');

channel_suffix  = ['_' channel '_tbytDat_dff'];
fileBeh         = GrabFiles_sort_trials(channel_suffix,0,{fullfile(filePath,'Matfiles')});
load(fullfile(fileBeh{1}),'tbytDat');

evt_file = GrabFiles_sort_trials('evtInS',0,GrabFiles_sort_trials('_g',0,{filePath}));
if isempty(evt_file)
    [evt_file,evt_folder] = uigetfile('*.mat','Select evtInS file',filePath);
    evt_file = fullfile(evt_folder,evt_file);
else
    evt_file = fullfile(evt_file{1});
end
load(evt_file,'evtInS');

% Re-insert manual water times into tbytDat
tbytDat = organizeManualWaterIntoTbytDat(evtInS.manualWater,tbytDat);

%% 2. Quick bookkeeping for manual vs. voluntary licks
nTrials     = numel(tbytDat);
mWaterV     = NaN(1,nTrials);
hitLickV    = NaN(1,nTrials);

hasMwater   = ~cellfun('isempty',{tbytDat.manualW});
mWaterV(hasMwater) = cellfun(@(a)a(1),{tbytDat(hasMwater).manualW});

hasHit      = ~cellfun('isempty',{tbytDat.hitLicks});
hitLickV(hasHit)   = cellfun(@(a)a(1),{tbytDat(hasHit).hitLicks});

voluntaryHitLickI  = (isnan(mWaterV) & ~isnan(hitLickV)) | (hitLickV < mWaterV);

%% 3. Load dF/F stacks and Allen masks
load(fullfile(filePath,'Matfiles',sprintf('%s_%s_dff_smCollect.mat',header,channel)), ...
     'dffsmCell');

load('allenDorsalMap.mat', ...
     'dorsalMaps','motorMask','smotorMask','ssMask','vMask','rsMask');

masks   = {motorMask, smotorMask, ssMask, vMask, rsMask};
regions = {'m1','m2','ss','v1','rs'};          % no hemi suffixes

% Transformation from imaging plane to AllenCCF
img_folder         = GrabFiles_sort_trials('_img',0,{filePath});
transParams_file   = find_keyword_containing_files(img_folder{1}, ...
                       ['transParamsAllen_' channel],'recursive',true);
load(fullfile(transParams_file{1}),'transParams');

%% 4. Pre-allocate result structure
eventFields  = {'stimOnDffC','postStimLickDff','firstWater', ...
                'hitLickFst','hitLickLast','hitLickToneAlign','consumeLick'};

rez = struct;
for ef = 1:numel(eventFields)
    for rl = 1:numel(regions)
        rez.(eventFields{ef}).(regions{rl}) = cell(nTrials,2);
    end
end

%% 5. Main processing loop
for t = 1:nTrials

    % --- Align current stack to Allen CCF -------------------------------
    dff_aligned = alignStackToAllenKabsch(dffsmCell{t}, ...
                                          dorsalMaps.dorsalMap, ...
                                          transParams.tformObj);

    % --- Extract region-wise mean traces --------------------------------
    dff = cell(1,numel(regions));
    for k = 1:numel(regions)
        dff{k} = apply2DMaskTo3DStack(dff_aligned, masks{k});
    end

    frT = tbytDat(t).frameTrel;          % frame time base

    % ---- Event alignments ---------------------------------------------
    rez = alignAndStoreEvent(dff,regions,rez,'stimOnDffC', ...
                             0,frT,[-0.9 5],t);

    % Post-stim lick bursts
    for jj = 1:numel(tbytDat(t).postStimChunk)
        evt_time = tbytDat(t).postStimChunk{jj}(1);
        if evt_time + 1 < frT(end)
            rez = alignAndStoreEvent(dff,regions,rez,'postStimLickDff', ...
                                     evt_time,frT,[-1 1],t,jj);
        end
    end

    % First water
    if ~isempty(tbytDat(t).water)
        firstWater = tbytDat(t).water(1) - tbytDat(t).evtOn;
        rez = alignAndStoreEvent(dff,regions,rez,'firstWater', ...
                                 firstWater,frT,[-1 1],t);
    end

    % Hit-aligned licks
    if ~isempty(tbytDat(t).hitLicks)
        rez = alignAndStoreEvent(dff,regions,rez,'hitLickFst', ...
                                 tbytDat(t).hitLicks(1),frT,[-1 1],t);
        rez = alignAndStoreEvent(dff,regions,rez,'hitLickLast', ...
                                 tbytDat(t).hitLicks(end),frT,[-1 1],t);
        rez = alignAndStoreEvent(dff,regions,rez,'hitLickToneAlign', ...
                                 0,frT,[-0.9 1],t);

        if ~isempty(tbytDat(t).consumeLicks)
            rez = alignAndStoreEvent(dff,regions,rez,'consumeLick', ...
                                     tbytDat(t).consumeLicks(1),frT,[-1 1],t);
        end
    end

    fprintf('processed trial #%d\n',t);
end

%% 6. Save
save(fullfile(filePath,'Matfiles', ...
     sprintf('%s_%s_dff_evtAligned_regionMask.mat',header,channel)), ...
     'rez','voluntaryHitLickI');
end
%%-----------------------------------------------------------------------
%% Helper functions
%-----------------------------------------------------------------------
function rez = alignAndStoreEvent(dffC,regions,rez, ...
                                  eventField,eventTime,frT,window, ...
                                  trialIdx,varargin)
% Align each region’s trace to a behavioural event and store it.
if ~isfield(rez,eventField)
    rez.(eventField) = struct();
end

for k = 1:numel(regions)
    regionName = regions{k};

    [alignedSig,alignedT] = alignToEvent(dffC{k},eventTime,frT,window);

    if ~isempty(varargin)                % e.g. postStim chunks
        jj = varargin{1};
        rez.(eventField).(regionName){trialIdx,1}{jj,1} = alignedSig;
        rez.(eventField).(regionName){trialIdx,1}{jj,2} = alignedT;
    else
        rez.(eventField).(regionName){trialIdx,1} = alignedSig;
        rez.(eventField).(regionName){trialIdx,2} = alignedT;
    end
end
end
%-----------------------------------------------------------------------
function [evtAlignedDffTs,evtAlignedTs] = ...
              alignToEvent(dffTs,eventTimeToAlign,frameT,timeWin)
% Extract a window around the event; return NaNs if window exceeds data.
win = timeWin + eventTimeToAlign;

frameI = frameT >= win(1) & frameT <= win(2);
evtAlignedDffTs = dffTs(frameI);
evtAlignedDffTs = evtAlignedDffTs(:).';          % row vector
evtAlignedTs    = frameT(frameI) - eventTimeToAlign;
evtAlignedTs    = evtAlignedTs(:).';
end
%-----------------------------------------------------------------------
function [meanTrace, imgStackMasked] = apply2DMaskTo3DStack(imgStack,mask2D)
% Mean dF/F within a 2-D mask across all frames.
mask3D          = mask2D(:,:,ones(1,size(imgStack,3)));  % replicate mask
imgStackMasked  = imgStack;
imgStackMasked(~mask3D) = NaN;
meanTrace       = squeeze(nanmean(nanmean(imgStackMasked,1),2));
end

