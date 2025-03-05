function trI = trialTypeInfoAuditoryGng(filePath)

match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

%% more indices for trial-type identification
trI.waterI = cellfun(@(a) ~isempty(a), {tbytDat.water}); 
trI.lickI = cellfun(@(a) ~isempty(a), {tbytDat.Lick}); 
trI.airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff}); 
trI.goI = [tbytDat.rewardTrI]'==1; 
trI.nogoI = [tbytDat.punishTrI]'==1; 
trI.hitI = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})'; 
trI.missI = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})'; 
trI.faI = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})'; 
trI.crI = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})'; 