function p = parse_input_ctxWideImagePreprocess(filePath, vargs)
    % parse input, and extract name-value pairs
    default_redoManual = false;         % logic to redo manual curation
    default_dffMethod = 'movingavg';    % default method for dff is moving average
    default_winF = 30;       % default window width (30 s) for F

    p = inputParser; % create parser object
    addRequired(p,'filePath')
    addParameter(p,'redoManual', default_redoManual)
    addParameter(p, 'dffMethod', default_dffMethod)
    addParameter(p, 'winF', default_winF)

    parse(p, filePath, vargs{:})
end