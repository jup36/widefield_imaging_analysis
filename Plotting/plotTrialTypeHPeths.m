
function plotTrialTypeHPeths(tbytDat_hAligned, trI, motifK)

timeX = getFirstNonEmptyTime(tbytDat_hAligned);

if isempty(timeX)
    warning('No valid H aligned time vector found.');
    return;
end

trialTypes = {};
trialMasks = {};

if isfield(trI, 'crI')
    trialTypes{end+1} = 'CR'; %#ok<AGROW>
    trialMasks{end+1} = trI.crI; %#ok<AGROW>
end

if isfield(trI, 'hitI')
    trialTypes{end+1} = 'Hit'; %#ok<AGROW>
    trialMasks{end+1} = trI.hitI; %#ok<AGROW>
end

if isfield(trI, 'missI')
    trialTypes{end+1} = 'Miss'; %#ok<AGROW>
    trialMasks{end+1} = trI.missI; %#ok<AGROW>
end

if isfield(trI, 'faI')
    trialTypes{end+1} = 'FA'; %#ok<AGROW>
    trialMasks{end+1} = trI.faI; %#ok<AGROW>
end

figure;
hold on;

for i = 1:numel(trialTypes)

    trialMask = logical(trialMasks{i});
    pethMat = extractMotifAlignedMatrix(tbytDat_hAligned, trialMask, motifK);

    if isempty(pethMat)
        continue;
    end

    plot(timeX, mean(pethMat, 1, 'omitnan'), 'LineWidth', 2);
end

xline(0, '--');
xlabel('Time from tone onset, sec');
ylabel(sprintf('Motif %d H', motifK));
title(sprintf('Trial-aligned H, motif %d', motifK));
legend(trialTypes, 'Location', 'best');
box off;

end