function plotMeanSemHweightings(trIdC, hCell, hTimeCell, tbytDat, blockIds, trainPulseEdgeC, motif_ind, timeWin, timeStep)

assert(isscalar(motif_ind)); % just work with h of one motif at a time for now 
tint = timeWin(1):timeStep:timeWin(2); 

hIntAlignedC = cell(1, length(trIdC)); 

for i = 1:length(trIdC)
    tempTrs = trIdC{i};
    hIntAlignedC{i} = cell(length(tempTrs), 1); 
    
    for t = 1:length(tempTrs) 
        hIntAlignedC{i}{t, 1} = NaN(1, length(tint)); 
        tr = tempTrs(t); 
        blockId = blockIds(tr); 
        framesOfTrial = cell2mat(trainPulseEdgeC{tr}); 
        framesOfTrial(end) = min(size(hCell{1, blockId}, 2), framesOfTrial(end));  
        hmatTr = hCell{2, blockId}(motif_ind, framesOfTrial(1):framesOfTrial(end));
        hTimeTr = hTimeCell{1, blockId}(1, framesOfTrial(1):framesOfTrial(end)) - tbytDat(tr).evtOn;
        tintI = hTimeTr(1) <= tint & tint <= hTimeTr(end); 
        
        hmatTrInt = interp1(hTimeTr, hmatTr, tint(tintI)); 
        hIntAlignedC{i}{t, 1}(1, tintI) = hmatTrInt; 
    end
end

% Convert nested cells to matrices
hIntAligned = cellfun(@cell2mat, hIntAlignedC, 'UniformOutput', false); 
mean_hIntAligned = cellfun(@nanmean, hIntAligned, 'UniformOutput', false); 
sem_hIntAligned = cellfun(@(x) nanstd(x) ./ sqrt(sum(~isnan(x), 1)), hIntAligned, 'UniformOutput', false); 

% Plot
figure; hold on;
colors = lines(length(mean_hIntAligned)); % Distinct colors for each group

for ii = 1:length(mean_hIntAligned)
    m = mean_hIntAligned{ii};
    s = sem_hIntAligned{ii};
    c = colors(ii, :);

    % Plot shaded error band
    fill([tint, fliplr(tint)], [m + s, fliplr(m - s)], c, ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Plot mean trace
    plot(tint, m, 'Color', c, 'LineWidth', 2);
end

xlabel('Time (s)');
ylabel('Signal');
title('Mean ± SEM of h-weightings');
legend(arrayfun(@(x) sprintf('Group %d', x), 1:length(trIdC), 'UniformOutput', false));
grid on;

end
