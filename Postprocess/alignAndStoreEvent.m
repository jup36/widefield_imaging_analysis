function rez = alignAndStoreEvent(dffC,regions,rez, ...
                                  eventField,eventTime,frT,window, ...
                                  trialIdx,varargin)
% Align each region’s trace to a behavioural event and store it.
if ~isfield(rez,eventField)
    rez.(eventField) = struct();
end

for k = 1:numel(regions)
    regionName = regions{k};

    [alignedSig,alignedT] = alignToEvent(dffC{k},eventTime,frT,window);

    if ~isempty(varargin)                % e.g. postStim chunks
        jj = varargin{1};
        rez.(eventField).(regionName){trialIdx,1}{jj,1} = alignedSig;
        rez.(eventField).(regionName){trialIdx,1}{jj,2} = alignedT;
    else
        rez.(eventField).(regionName){trialIdx,1} = alignedSig;
        rez.(eventField).(regionName){trialIdx,2} = alignedT;
    end
end
end