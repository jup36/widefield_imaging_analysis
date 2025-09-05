function hIntAlignedC = tbytHweightingsPsth(hCell, tbytDat, timeWin, timeStep)
% timeWin = [-0.9 5];
% timeStep = 0.01;

blockIds = [tbytDat.limeLEDTrainI]';
trainPulseEdgeC = {tbytDat.limeLEDPulsesOfTrain};

tint = timeWin(1):timeStep:timeWin(2);

hIntAlignedC = cell(2, numel(tbytDat));

for t = 1:numel(tbytDat)
    blockId = blockIds(t);
    if blockId <= size(hCell, 2)
        framesOfTrial = cell2mat(trainPulseEdgeC{t});
        framesOfTrial(end) = min(size(hCell{1, blockId}, 2), framesOfTrial(end));
        hmatTr = hCell{2, blockId}(:, framesOfTrial(1):framesOfTrial(end));
        hTimeTr = hCell{3, blockId}(1, framesOfTrial(1):framesOfTrial(end)) - tbytDat(t).evtOn;
        tintI = hTimeTr(1) <= tint & tint <= hTimeTr(end);

        % Preallocate
        hmatTrInt = zeros(size(hmatTr, 1), sum(tintI));

        % Interpolate each motif (row) over the new time points
        for j = 1:size(hmatTr, 1)
            hmatTrInt(j, :) = interp1(hTimeTr, hmatTr(j, :), tint(tintI), 'linear', 'extrap');
        end
        hIntAlignedC{1, t} = hmatTrInt;
        hIntAlignedC{2, t} = tint(tintI);
    end
end
