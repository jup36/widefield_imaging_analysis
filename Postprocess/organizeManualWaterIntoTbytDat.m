function tbytDat = organizeManualWaterIntoTbytDat(mwts, tbytDat)

% get relative manual water time (relative to the tone onset)
relTimeC = cellfun(@(a) mwts-a, {tbytDat.evtOn}, 'UniformOutput', false); 
% logic indicating the trial membership of each manual water reward
manualWaterLogicC = cellfun(@(a) mwts-a>0 & mwts-a<6, {tbytDat.evtOn}, 'UniformOutput', false);
% relative manual water timestamps per trial
relTimeManualWaterC = cellfun(@(a, b) a(b), relTimeC, manualWaterLogicC, 'UniformOutput', false); 
% add the info to tbytDat
[tbytDat.manualW] = relTimeManualWaterC{:};

end