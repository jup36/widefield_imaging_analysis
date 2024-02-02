function p = parse_input_ctxWideImagePreprocess(filePath, vargs)
    % parse input, and extract name-value pairs
    default_redoManual = false;    % logic to redo manual curation

    p = inputParser; % create parser object
    addRequired(p, 'filePath')
    addParameter(p, 'redoManual', default_redoManual)

    parse(p, filePath, vargs{:})
end
