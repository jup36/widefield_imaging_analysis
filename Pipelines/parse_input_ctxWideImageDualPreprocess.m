    function p = parse_input_ctxWideImageDualPreprocess(filePath, vargs)
        % parse input, and extract name-value pairs
        default_redoManual = false;         % logic to redo manual curation
        default_dffMethod = 'movingavg';    % default method for dff is moving average
        default_winF = 30;                  % default window width (30 s) for F
        default_channelToProcess = 'both';
        default_byPassManualAlign = false;  % if true, skip ManualAlignmentAdjust and redo vasculature only

        p = inputParser;
        addRequired(p,'filePath')
        addParameter(p,'redoManual', default_redoManual)
        addParameter(p, 'dffMethod', default_dffMethod)
        addParameter(p, 'winF', default_winF)
        addParameter(p, 'channelToProcess', default_channelToProcess)
        addParameter(p, 'byPassManualAlign', default_byPassManualAlign)

        parse(p, filePath, vargs{:})
    end
