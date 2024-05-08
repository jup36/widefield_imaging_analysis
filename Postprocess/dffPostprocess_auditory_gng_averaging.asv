function dffPostprocess_auditory_gng_averaging(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';

%% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh)
    fileBeh = GrabFiles_sort_trials('tbytDat_dff', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

% get behavioral data further analyzed
fileBehParseGng = GrabFiles_sort_trials('tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehParseGng)
    tbytDat = parseAuditoryGngTrials(tbytDat);
    save(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')
else
    load(fileBehParseGng{1})
end

% load the preprocessed dffs collect
load(fullfile(filePath, 'Matfiles', strcat(header, '_dffsmCollect.mat')), 'dffsmCell');

% load Allen dorsalMap
load(fullfile('/Users/jp3025/Documents/codes/WidefieldAnalysis_Musall/allenDorsalMap.mat'), 'dorsalMaps', 'motorMask', 'smotorMask', 'ssMask', 'vMask', 'rsMask');

% load Allen transformation object
folder_list_imgTrial = GrabFiles_sort_trials(header, 0, {fullfile(filePath, strcat(header, '_img'))}); % use GrabFiles_sort_trials to sort both files and folders
load(fullfile(folder_list_imgTrial{1}, 'transParamsAllen.mat'), 'transParams')

% file directory for trials
filePathTrials = fullfile(filePath, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

for t = 1:length(tbytDat)

    % align image stack to the AllenCCF
    dff = alignStackToAllenKabsch(dffsmCell{t}, dorsalMaps.dorsalMap, transParams.tformObj); % row x col x frame

    frT = tbytDat(t).frameTrel; % store timestamps

    % apply region masks
    dffM1 = apply2DMaskTo3DStack(dff, motorMask);
    dffM2 = apply2DMaskTo3DStack(dff, smotorMask);
    dffSs = apply2DMaskTo3DStack(dff, ssMask);
    dffV1 = apply2DMaskTo3DStack(dff, vMask);
    dffRs = apply2DMaskTo3DStack(dff, rsMask);

    %% Align to tone Onset
    [rez.stimOnDffC.m1{t, 1}, rez.stimOnDffC.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-0.9 5]);
    [rez.stimOnDffC.m2{t, 1}, rez.stimOnDffC.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-0.9 5]);
    [rez.stimOnDffC.ss{t, 1}, rez.stimOnDffC.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-0.9 5]);
    [rez.stimOnDffC.v1{t, 1}, rez.stimOnDffC.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-0.9 5]);
    [rez.stimOnDffC.rs{t, 1}, rez.stimOnDffC.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-0.9 5]);

    %% Align to licks during ITI (technically impossible with the retractable spout)
    %     for ii = 1:length(tbytDat(t).itiLickChunk) % there can be multiple bouts
    %         [rez.itiLickDff.m1{t, 1}{ii, 1}, rez.itiLickDff.m1{t, 1}{ii, 2}] = alignToEvent(dffM1, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % M1 align to the ITI lick bout use the 1st lick of each bout
    %         [rez.itiLickDff.m2{t, 1}{ii, 1}, rez.itiLickDff.m2{t, 1}{ii, 2}]  = alignToEvent(dffM2, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % M2 align to the ITI lick bout use the 1st lick of each bout
    %         [rez.itiLickDff.ss{t, 1}{ii, 1}, rez.itiLickDff.ss{t, 1}{ii, 2}] = alignToEvent(dffSs, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % Ss align to the ITI lick bout use the 1st lick of each bout
    %         [rez.itiLickDff.v1{t, 1}{ii, 1}, rez.itiLickDff.v1{t, 1}{ii, 2}] = alignToEvent(dffV1, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % V1 align to the ITI lick bout use the 1st lick of each bout
    %         [rez.itiLickDff.rs{t, 1}{ii, 1}, rez.itiLickDff.rs{t, 1}{ii, 2}] = alignToEvent(dffRs, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % Rs align to the ITI lick bout use the 1st lick of each bout
    %     end

    %% Align to postStim licks after tone offset
    for jj = 1:length(tbytDat(t).postStimChunk)
        if tbytDat(t).postStimChunk{jj}(1) + 1 < frT(end)
            [rez.postStimLickDff.m1{t, 1}{jj, 1}, rez.postStimLickDff.m1{t, 1}{jj, 2}] = alignToEvent(dffM1, tbytDat(t).postStimChunk{jj}(1), frT, [-1 1]); % M1 align to the postStimLick lick bout use the 1st lick of each bout
            [rez.postStimLickDff.m2{t, 1}{jj, 1}, rez.postStimLickDff.m2{t, 1}{jj, 2}] = alignToEvent(dffM2, tbytDat(t).postStimChunk{jj}(1), frT, [-1 1]); % M2 align to the postStimLick lick bout use the 1st lick of each bout
            [rez.postStimLickDff.ss{t, 1}{jj, 1}, rez.postStimLickDff.ss{t, 1}{jj, 2}] = alignToEvent(dffSs, tbytDat(t).postStimChunk{jj}(1), frT, [-1 1]); % Ss align to the postStimLick lick bout use the 1st lick of each bout
            [rez.postStimLickDff.v1{t, 1}{jj, 1}, rez.postStimLickDff.v1{t, 1}{jj, 2}] = alignToEvent(dffV1, tbytDat(t).postStimChunk{jj}(1), frT, [-1 1]); % V1 align to the postStimLick lick bout use the 1st lick of each bout
            [rez.postStimLickDff.rs{t, 1}{jj, 1}, rez.postStimLickDff.rs{t, 1}{jj, 2}] = alignToEvent(dffRs, tbytDat(t).postStimChunk{jj}(1), frT, [-1 1]); % Rs align to the postStimLick lick bout use the 1st lick of each bout
        end
    end

    %% xcorr dff and lick bouts
    % lick chuncks time relative to stimOn
    if ~isempty(tbytDat(t).LickChunk)
        tbytDat(t).LickBoutRel = cellfun(@(a) a-tbytDat(t).evtOn, tbytDat(t).LickChunk, 'UniformOutput', false);
        %m1
        [lickOnTf, lickBoxOnTf, dffsOnTfItpM1, ~] = alignDffAndLickBouts(tbytDat(t).LickBoutRel, dffM1, tbytDat(1).frameTrel, [-0.9 5]);
        rez.lickOnTfBin{t, 1} = bin1msSpkCountMat( lickOnTf, 50, 50 );
        rez.dffsOnTfItpM1{t, 1} = dffsOnTfItpM1(1:50:end);
        [rez.xcorrDffLick.m1{t, 1}, rez.xcorrDffLick.m1{t, 2}] = xcorr(dffsOnTfItpM1-min(dffsOnTfItpM1), lickOnTf, 2000, 'normalized');
        [rez.xcorrDffLickBox.m1{t, 1}, rez.xcorrDffLickBox.m1{t, 2}] = xcorr(dffsOnTfItpM1-min(dffsOnTfItpM1), lickBoxOnTf, 2000, 'normalized');
        [rez.xcorrDsDffLickBox.m1{t, 1}, rez.xcorrDsDffLickBox.m1{t, 2}] = xcorr(dffsOnTfItpM1(1:50:end)-min(dffsOnTfItpM1), lickBoxOnTf(1:50:end), 40, 'normalized');
        %m2
        [lickOnTf, lickBoxOnTf, dffsOnTfItpM2, ~] = alignDffAndLickBouts(tbytDat(t).LickBoutRel, dffM2, tbytDat(1).frameTrel, [-0.9 5]);
        rez.dffsOnTfItpM2{t, 1} = dffsOnTfItpM2(1:50:end);
        [rez.xcorrDffLick.m2{t, 1}, rez.xcorrDffLick.m2{t, 2}] = xcorr(dffsOnTfItpM2-min(dffsOnTfItpM2), lickOnTf, 2000, 'normalized');
        [rez.xcorrDffLickBox.m2{t, 1}, rez.xcorrDffLickBox.m2{t, 2}] = xcorr(dffsOnTfItpM2-min(dffsOnTfItpM2), lickBoxOnTf, 2000, 'normalized');
        [rez.xcorrDsDffLickBox.m2{t, 1}, rez.xcorrDsDffLickBox.m2{t, 2}] = xcorr(dffsOnTfItpM2(1:50:end)-min(dffsOnTfItpM2), lickBoxOnTf(1:50:end), 2000, 'normalized');
        %ss
        [lickOnTf, lickBoxOnTf, dffsOnTfItpSs, ~] = alignDffAndLickBouts(tbytDat(t).LickBoutRel, dffSs, tbytDat(1).frameTrel, [-0.9 5]);
        rez.dffsOnTfItpSs{t, 1} = dffsOnTfItpSs(1:50:end);
        [rez.xcorrDffLick.ss{t, 1}, rez.xcorrDffLick.ss{t, 2}] = xcorr(dffsOnTfItpSs-min(dffsOnTfItpSs), lickOnTf, 2000, 'normalized');
        [rez.xcorrDffLickBox.ss{t, 1}, rez.xcorrDffLickBox.ss{t, 2}] = xcorr(dffsOnTfItpSs-min(dffsOnTfItpSs), lickBoxOnTf, 2000, 'normalized');
        [rez.xcorrDsDffLickBox.ss{t, 1}, rez.xcorrDsDffLickBox.ss{t, 2}] = xcorr(dffsOnTfItpSs(1:50:end)-min(dffsOnTfItpSs), lickBoxOnTf(1:50:end), 2000, 'normalized');
        %v1
        [lickOnTf, lickBoxOnTf, dffsOnTfItpV1, ~] = alignDffAndLickBouts(tbytDat(t).LickBoutRel, dffV1, tbytDat(1).frameTrel, [-0.9 5]);
        rez.dffsOnTfItpV1{t, 1} = dffsOnTfItpV1(1:50:end);
        [rez.xcorrDffLick.V1{t, 1}, rez.xcorrDffLick.V1{t, 2}] = xcorr(dffsOnTfItpV1-min(dffsOnTfItpV1), lickOnTf, 2000, 'normalized');
        [rez.xcorrDffLickBox.v1{t, 1}, rez.xcorrDffLickBox.v1{t, 2}] = xcorr(dffsOnTfItpV1-min(dffsOnTfItpV1), lickBoxOnTf, 2000, 'normalized');
        [rez.xcorrDsDffLickBox.v1{t, 1}, rez.xcorrDsDffLickBox.v1{t, 2}] = xcorr(dffsOnTfItpV1(1:50:end)-min(dffsOnTfItpV1), lickBoxOnTf(1:50:end), 2000, 'normalized');
        %rs
        [lickOnTf, lickBoxOnTf, dffsOnTfItpRs, ~] = alignDffAndLickBouts(tbytDat(t).LickBoutRel, dffRs, tbytDat(1).frameTrel, [-0.9 5]);
        rez.dffsOnTfItpRs{t, 1} = dffsOnTfItpRs(1:50:end);
        [rez.xcorrDffLick.rs{t, 1}, rez.xcorrDffLick.rs{t, 2}] = xcorr(dffsOnTfItpRs-min(dffsOnTfItpRs), lickOnTf, 2000, 'normalized');
        [rez.xcorrDffLickBox.rs{t, 1}, rez.xcorrDffLickBox.rs{t, 2}] = xcorr(dffsOnTfItpRs-min(dffsOnTfItpRs), lickBoxOnTf, 2000, 'normalized');
        [rez.xcorrDsDffLickBox.rs{t, 1}, rez.xcorrDsDffLickBox.rs{t, 2}] = xcorr(dffsOnTfItpRs(1:50:end)-min(dffsOnTfItpRs), lickBoxOnTf(1:50:end), 2000, 'normalized');
    end
    %figure; hold on;
    %plot(tF, smooth2a(dffsOnTfItp, 0, 50));
    %plot(tF, lickBoxOnTf);
    %plot(tF, licksOnTf);
    %set(gca, 'XTick', -5:1:5, 'TickDir', 'out', 'YTick', -1.5:0.5:1.5);
    %title('Dff overlayed with licks')
    %print(fullfile(filePath, 'Figure', 'dff_overlayed_with_licks.pdf'), '-dpdf', '-vector', '-bestfit')

    %% Align to water (1st water only)
    if ~isempty(tbytDat(t).water)
        firstWater = tbytDat(t).water(1) - tbytDat(t).evtOn;
        % m1
        [rez.firstWater.m1{t, 1}, rez.firstWater.m1{t, 1}] = alignToEvent(dffM1, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % m2
        [rez.firstWater.m2{t, 1}, rez.firstWater.m2{t, 2}] = alignToEvent(dffM2, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % ss
        [rez.firstWater.ss{t, 1}, rez.firstWater.ss{t, 2}] = alignToEvent(dffSs, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % v1
        [rez.firstWater.v1{t, 1}, rez.firstWater.v1{t, 2}] = alignToEvent(dffV1, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % rs
        [rez.firstWater.rs{t, 1}, rez.firstWater.rs{t, 2}]  = alignToEvent(dffRs, firstWater, frT, [-1 1]); % consume lick 1st in the bout
    end

    %% Align to airpuff (1st airpuff only)
    if ~isempty(tbytDat(t).airpuff)
        firstAirpuff = tbytDat(t).airpuff(1) - tbytDat(t).evtOn;
        % m1
        [rez.firstAirpuff.m1{t, 1}, rez.firstAirpuff.m1{t, 2}] = alignToEvent(dffM1, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % m2
        [rez.firstAirpuff.m2{t, 1}, rez.firstAirpuff.m2{t, 2}] = alignToEvent(dffM2, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % ss
        [rez.firstAirpuff.ss{t, 1}, rez.firstAirpuff.ss{t, 2}] = alignToEvent(dffSs, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % v1
        [rez.firstAirpuff.v1{t, 1}, rez.firstAirpuff.v1{t, 2}] = alignToEvent(dffV1, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % rs
        [rez.firstAirpuff.rs{t, 1}, rez.firstAirpuff.rs{t, 2}]  = alignToEvent(dffRs, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
    end

    %% Align to licks during Cue
    % Hit
    % ToDo: Align to the first or last lick?
    if ~isempty(tbytDat(t).hitLicks)
        % m1
        [rez.hitDffFirstLick.m1{t, 1}, rez.hitDffFirstLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).hitLicks(1), frT, [-1 1]); % hit lick 1st in the bout
        [rez.hitDffLastLick.m1{t, 1}, rez.hitDffLastLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).hitLicks(end), frT, [-1 1]); % hit lick last in the bout
        [rez.hitDffCueOn.m1{t, 1}, rez.hitDffCueOn.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-0.9 1]); % hit lick cueOn
        % m2
        [rez.hitDffFirstLick.m2{t, 1}, rez.hitDffFirstLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.m2{t, 1}, rez.hitDffLastLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.m2{t, 1}, rez.hitDffCueOn.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-0.9 1]);
        % ss
        [rez.hitDffFirstLick.ss{t, 1}, rez.hitDffFirstLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.ss{t, 1}, rez.hitDffLastLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.ss{t, 1}, rez.hitDffCueOn.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-0.9 1]);
        % v1
        [rez.hitDffFirstLick.v1{t, 1}, rez.hitDffFirstLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.v1{t, 1}, rez.hitDffLastLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.v1{t, 1}, rez.hitDffCueOn.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-0.9 1]);
        % rs
        [rez.hitDffFirstLick.rs{t, 1}, rez.hitDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.rs{t, 1}, rez.hitDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.rs{t, 1}, rez.hitDffCueOn.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-0.9 1]);

        if ~isempty(tbytDat(t).consumeLicks)
            % m1
            [rez.hitDffConsumeLick.m1{t, 1}, rez.hitDffConsumeLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % m2
            [rez.hitDffConsumeLick.m2{t, 1}, rez.hitDffConsumeLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % ss
            [rez.hitDffConsumeLick.ss{t, 1}, rez.hitDffConsumeLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % v1
            [rez.hitDffConsumeLick.v1{t, 1}, rez.hitDffConsumeLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % rs
            [rez.hitDffConsumeLick.rs{t, 1}, rez.hitDffConsumeLick.rs{t, 2}]  = alignToEvent(dffRs, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
        end

    % Miss
    elseif tbytDat(t).rewardTrI && isempty(tbytDat(t).water)
        % m1
        [rez.missCueDff.m1{t, 1}, rez.missCueDff.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-0.9 1]); % align to cue onset
        % m2
        [rez.missCueDff.m2{t, 1}, rez.missCueDff.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-0.9 1]); % align to cue onset
        % ss
        [rez.missCueDff.ss{t, 1}, rez.missCueDff.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-0.9 1]); % align to cue onset
        % v1
        [rez.missCueDff.v1{t, 1}, rez.missCueDff.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-0.9 1]); % align to cue onset
        % rs
        [rez.missCueDff.rs{t, 1}, rez.missCueDff.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-0.9 1]); % align to cue onset

        % False Alarm
    elseif ~isempty(tbytDat(t).faLicks)
        % m1
        [rez.faDffFirstLick.m1{t, 1}, rez.faDffFirstLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).faLicks(1), frT, [-1 1]); % fa lick 1st in the bout
        [rez.faDffLastLick.m1{t, 1}, rez.faDffLastLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).faLicks(end), frT, [-1 1]); % fa lick last in the bout
        [rez.faDffCueOn.m1{t, 1}, rez.faDffCueOn.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-0.9 1]); % fa lick cueOn
        % m2
        [rez.faDffFirstLick.m2{t, 1}, rez.faDffFirstLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.m2{t, 1}, rez.faDffLastLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.m2{t, 1}, rez.faDffCueOn.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-0.9 1]);
        % ss
        [rez.faDffFirstLick.ss{t, 1}, rez.faDffFirstLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.ss{t, 1}, rez.faDffLastLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.ss{t, 1}, rez.faDffCueOn.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-0.9 1]);
        % v1
        [rez.faDffFirstLick.v1{t, 1}, rez.faDffFirstLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.v1{t, 1}, rez.faDffLastLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.v1{t, 1}, rez.faDffCueOn.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-0.9 1]);
        % rs
        [rez.faDffFirstLick.rs{t, 1}, rez.faDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.rs{t, 1}, rez.faDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.rs{t, 1}, rez.faDffCueOn.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-0.9 1]);

        if ~isempty(tbytDat(t).postAirLicks)
            % m1
            [rez.faDffPostAirLick.m1{t, 1}, rez.faDffPostAirLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % m2
            [rez.faDffPostAirLick.m2{t, 1}, rez.faDffPostAirLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % ss
            [rez.faDffPostAirLick.ss{t, 1}, rez.faDffPostAirLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % v1
            [rez.faDffPostAirLick.v1{t, 1}, rez.faDffPostAirLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % rs
            [rez.faDffPostAirLick.rs{t, 1}, rez.faDffPostAirLick.rs{t, 2}]  = alignToEvent(dffRs, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
        end

        % Correct Rejection
    elseif tbytDat(t).punishTrI && isempty(tbytDat(t).airpuff)
        % m1
        [rez.crCueDff.m1{t, 1}, rez.crCueDff.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-0.9 1]); % align to cue onset
        % m2
        [rez.crCueDff.m2{t, 1}, rez.crCueDff.m2{t, 2}]  = alignToEvent(dffM2, 0, frT, [-0.9 1]); % align to cue onset
        % ss
        [rez.crCueDff.ss{t, 1}, rez.crCueDff.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-0.9 1]); % align to cue onset
        % v1
        [rez.crCueDff.v1{t, 1}, rez.crCueDff.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-0.9 1]); % align to cue onset
        % rs
        [rez.crCueDff.rs{t, 1}, rez.crCueDff.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-0.9 1]); % align to cue onset
    end
    fprintf("processed trial #%d\n", t)
end

% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [evtAlignedDffTs, evtAlignedTs] = alignToEvent(dffTs, eventTimeToAlign, frameT, timeWin)
        win = timeWin + eventTimeToAlign;

        if min(win) >= min(frameT) && max(win) <= max(frameT)
            frameI = frameT >= min(win) & frameT <= max(win);
            evtAlignedDffTs = dffTs(frameI);
            evtAlignedDffTs = evtAlignedDffTs(:).'; % make it a row vector
            evtAlignedTs = frameT(frameI);
            evtAlignedTs = evtAlignedTs(:).'; % make it a row vector
        else
            warning("The peri-event window goes out of bound!")
            evtAlignedDffTs = [];
            evtAlignedTs = [];
        end
    end


    function [imagStackMaskedMean, imgStackMasked] = apply2DMaskTo3DStack(imgStack, mask2D)
        % Replicate the 2D mask to match the 3D stack dimensions
        replicatedMask = repmat(mask2D, [1, 1, size(imgStack, 3)]);

        % Apply the mask: replace values in the stack where the mask is zero (or false) with NaN
        imgStack(~replicatedMask) = NaN;

        imgStackMasked = imgStack;

        imagStackMaskedMean = squeeze(nanmean(nanmean(imgStackMasked, 1), 2));

    end


    function [mRez, sRez] = trialGroupMeanSem(tbytPsthC, groupC)
        % 'tbytPsthC' contains trial-by-trial data in each cell (it is assumed that data are temporally aligned across trials already).
        % 'groupC' contains logical for each grouping, logicals are assumed to be of same lengths as the number of trials.
        % 'mRez' returns the group mean for each group logical per row.
        % 'sRez' returns the group sem for each group logical per row.

        % tbytPsthC = stimOnDff_m1_itp;
        % groupC = {[tbytDat.rewardTrI], [tbytDat.punishTrI]};

        % sanity check 1: all trial numbers must match!
        lenT = length(tbytPsthC);
        trialIC = cell2mat(cellfun(@length, groupC, 'UniformOutput', false));
        if length(unique([lenT, trialIC]))~=1
            error("The length of trials and length(s) of group logics do not match!")
        end

        % sanity check 2:
        if length(unique(cell2mat(cellfun(@length, tbytPsthC, 'UniformOutput', false))))~=1
            error("Some trials have different lengths, e.g., temporally misaligned!")
        end

        nGroups = length(groupC);

        mRez = [];
        sRez = [];

        for gg = 1:nGroups
            gDat = cell2mat(tbytPsthC(logical(groupC{gg}), 1));
            [gMean, ~, gSem] = meanstdsem(gDat);
            mRez = [mRez; gMean];
            sRez = [sRez; gSem];
        end
    end


    function accuracy = logireg(X, Y)
        % % X: observations-by-feature
        % %   e.g., cell2mat(stimOnDff_m1_itp);
        % % Y: label

        Y = Y(:);

        % Determine the number of samples for each class in the test set
        numClass0 = sum(Y == 0);
        numClass1 = sum(Y == 1);
        testSetSize = min(numClass0, numClass1) / 2; % Or any other criterion

        % Randomly sample for the test set
        indices0 = find(Y == 0);
        indices1 = find(Y == 1);
        testIndices0 = randsample(indices0, testSetSize);
        testIndices1 = randsample(indices1, testSetSize);

        % Combine test samples and create the test set
        testInd = [testIndices0; testIndices1];
        testInd = testInd(randperm(length(testInd)));

        % Use the remaining data as the training set
        trainInd = setdiff(1:length(Y), testInd);

        accuracy = nan(1, size(X, 2));
        % Training
        for ff = 1:size(X, 2) % iterate features
            model = fitglm(X(trainInd, ff), Y(trainInd), 'Distribution', 'binomial');
            % Testing
            Y_pred = predict(model, X(testInd,ff)) > 0.5; % Thresholding at 0.5
            accuracy(1, ff) = sum(Y_pred == Y(testInd)) / length(Y_pred);
        end

    end

    function [licksOnTf, lickBoxOnTf, dffsOnTfItp, tF] = alignDffAndLickBouts(lickBouts, dffs, frameT, timeW)
        % lickBouts = tbytDat(1).LickChunkRel;
        % dffs = dffM1;
        % frameT = tbytDat(1).frameTrel;
        % timeW = [-5 5];

        % get licks on the specified time frame
        tF = timeW(1):0.001:timeW(end); % time frame
        licksOnTf = zeros(1, length(tF));
        lickBoxOnTf = zeros(1, length(tF));

        for ts = 1:length(lickBouts)
            lb = lickBouts{ts};

            vlicks = lb(cell2mat(arrayfun(@(a) a>=min(tF) & a<=max(tF), lb, 'UniformOutput', false)));

            if ~isempty(vlicks)
                licksOnTfI = locateOnArray(tF, vlicks);
                licksOnTf(licksOnTfI) = 1;
                lickBoxOnTf(min(licksOnTfI):max(licksOnTfI)) = 1;
            end
        end

        % get dff on the specified time frame with interpolation
        ftI = frameT >= min(tF) & frameT <= max(tF);
        dffsOnTf = dffs(ftI);
        frameTOnTf = frameT(ftI);
        % Interpolate
        dffsOnTfItp = interp1(frameTOnTf, dffsOnTf, tF, 'linear', 'extrap');

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function loc = locateOnArray(array, target)
            % Convert the array to a column vector for broadcasting
            array = array(:);
            target = target(:);

            % Subtract each target value from the entire array and find the minimum absolute value
            [~, loc] = min(abs(array - target.'), [], 1);
        end


    end



end


