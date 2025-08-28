function trI = trialTypeInfoAuditoryGng(filePath)
filePathC = GrabFiles_sort_trials('_tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
if isempty(filePathC)
    error("Behavioral data (tbytDat_parseGng) not available!")
else
    load(filePathC{1}, 'tbytDat')

    %% more indices for trial-type identification
    trI.waterI = cellfun(@(a) ~isempty(a), {tbytDat.water});
    trI.lickI = cellfun(@(a) ~isempty(a), {tbytDat.Lick});
    trI.airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff});
    trI.goI = [tbytDat.rewardTrI]'==1;
    trI.nogoI = [tbytDat.punishTrI]'==1;
    trI.hitI = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})' & trI.waterI';
    trI.missI = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})';
    trI.crI = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})';
    if isfield(tbytDat, 'faLicks')
        trI.faI = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})' & trI.airpuffI';
    else
        trI.faI = cell2mat({tbytDat.punishTrI})' & trI.airpuffI';
    end
end