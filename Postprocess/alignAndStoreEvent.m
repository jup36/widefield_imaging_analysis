function rez = alignAndStoreEvent(dffL, dffR, regions, rez, eventField, eventTime, frT, window, trialIdx, varargin)

    if ~isfield(rez, eventField)
        rez.(eventField) = struct();
    end

    for k = 1:numel(regions)
        regionL = [regions{k}, 'L'];
        regionR = [regions{k}, 'R'];

        [alignedL, alignedTL] = alignToEvent(dffL{k}, eventTime, frT, window);
        [alignedR, alignedTR] = alignToEvent(dffR{k}, eventTime, frT, window);

        if ~isempty(varargin)
            jj = varargin{1};
            rez.(eventField).(regionL){trialIdx, 1}{jj, 1} = alignedL;
            rez.(eventField).(regionL){trialIdx, 1}{jj, 2} = alignedTL;
            rez.(eventField).(regionR){trialIdx, 1}{jj, 1} = alignedR;
            rez.(eventField).(regionR){trialIdx, 1}{jj, 2} = alignedTR;
        else
            rez.(eventField).(regionL){trialIdx, 1} = alignedL;
            rez.(eventField).(regionL){trialIdx, 2} = alignedTL;
            rez.(eventField).(regionR){trialIdx, 1} = alignedR;
            rez.(eventField).(regionR){trialIdx, 2} = alignedTR;
        end
    end
end
