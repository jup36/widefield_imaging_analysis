function pevIntAlignedC = tbytPerFramePevPsth(PevC, tbytDat, timeWin, timeStep)
% timeWin = [-0.9 5];
% timeStep = 0.01;

blockIds = [tbytDat.limeLEDTrainI]';
trainPulseEdgeC = {tbytDat.limeLEDPulsesOfTrain};

tint = timeWin(1):timeStep:timeWin(2);

pevIntAlignedC = cell(2, numel(tbytDat));

for t = 1:numel(tbytDat)
    blockId = blockIds(t);
    if blockId <= size(PevC, 2)
        framesOfTrial = cell2mat(trainPulseEdgeC{t});
        framesOfTrial(end) = min(size(PevC{1, blockId}, 2), framesOfTrial(end));
        PevmatTr = PevC{1, blockId}(:, framesOfTrial(1):framesOfTrial(end));
        PevTimeTr = PevC{2, blockId}(1, framesOfTrial(1):framesOfTrial(end)) - tbytDat(t).evtOn;
        %tintI = PevTimeTr(1) <= tint & tint <= PevTimeTr(end);

        % % Preallocate
        % PevmatTrInt = zeros(size(PevmatTr, 1), sum(tintI));
        % 
        % % Interpolate each motif (row) over the new time points
        % for j = 1:size(PevmatTr, 1)
        %     PevmatTrInt(j, :) = interp1(PevTimeTr, PevmatTr(j, :), tint(tintI), 'linear', 'extrap');
        % end
        % pevIntAlignedC{1, t} = PevmatTrInt;
        % pevIntAlignedC{2, t} = tint(tintI);

        PevmatTrInt = zeros(size(PevmatTr, 1), numel(tint));  % full length

        for j = 1:size(PevmatTr, 1)
            PevmatTrInt(j, :) = interp1(PevTimeTr, PevmatTr(j, :), tint, 'linear', 'extrap');
        end

        pevIntAlignedC{1, t} = PevmatTrInt;
        pevIntAlignedC{2, t} = tint;

    end
end
end