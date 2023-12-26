%% mean classification accuracy (hit, fa, iti)
rez.meanLickHitFaRestNb.m1 = mean(rez.lickHitFaRestNb.m1); 
rez.meanLickHitFaRestNb.m2 = mean(rez.lickHitFaRestNb.m2); 
rez.meanLickHitFaRestNb.ss = mean(rez.lickHitFaRestNb.ss); 
rez.meanLickHitFaRestNb.rs = mean(rez.lickHitFaRestNb.rs);
rez.meanLickHitFaRestNb.v1 = mean(rez.lickHitFaRestNb.v1);


x = -1:0.05:1; 

plotMeanSem(smooth2a([rez.meanLickHitFaRestNb.m1; rez.meanLickHitFaRestNb.m2; rez.meanLickHitFaRestNb.ss; rez.meanLickHitFaRestNb.rs; rez.meanLickHitFaRestNb.v1], 0, 1), ...
    zeros(5, size(rez.meanLickHitFaRestNb.m1, 2)), x, {'M1', 'M2', 'S1', 'Rs', 'V1'});
title("naive bayes classifier (Hit, FA, ITI)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.05:3); xlim([-1 1]); axis tight; grid on; 
print(fullfile(filePath, 'Figure', 'lick_classifier_naive_bayes_Hit_FA_ITI_across_regions.pdf'), '-dpdf', '-vector', '-bestfit')

%% mean classification accuracy (hit, fa)
rez.meanLickHitFaNb.m1 = mean(rez.lickHitFaNb.m1); 
rez.meanLickHitFaNb.m2 = mean(rez.lickHitFaNb.m2); 
rez.meanLickHitFaNb.ss = mean(rez.lickHitFaNb.ss); 
rez.meanLickHitFaNb.rs = mean(rez.lickHitFaNb.rs);
rez.meanLickHitFaNb.v1 = mean(rez.lickHitFaNb.v1);

plotMeanSem(smooth2a([rez.meanLickHitFaNb.m1; rez.meanLickHitFaNb.m2; rez.meanLickHitFaNb.ss; rez.meanLickHitFaNb.rs; rez.meanLickHitFaNb.v1], 0, 1), ...
    zeros(5, size(rez.meanLickHitFaNb.m1, 2)), x, {'M1', 'M2', 'S1', 'Rs', 'V1'});
title("naive bayes classifier (Hit, FA)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.05:3); xlim([-1 1]); axis tight; grid on; 
print(fullfile(filePath, 'Figure', 'lick_classifier_naive_bayes_Hit_FA_across_regions.pdf'), '-dpdf', '-vector', '-bestfit')