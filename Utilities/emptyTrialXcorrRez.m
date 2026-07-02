function rez = emptyTrialXcorrRez(maxLagBins, dt)

lagBins = -maxLagBins:maxLagBins;

rez = struct();
rez.rTrial = [];
rez.meanR = [];
rez.semR = [];
rez.peakR = [];
rez.peakLagSec = [];
rez.lagBins = lagBins;
rez.lagSec = lagBins * dt;
rez.validTrials = [];
rez.nTrials = 0;
rez.nTrialsPerMotif = [];
rez.subtractConditionMean = false;

end
