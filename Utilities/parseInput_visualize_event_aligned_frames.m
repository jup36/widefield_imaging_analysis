    function p = parseInput_visualize_event_aligned_frames(filePath, numPeriFrames, trial, vargs)
        % parse input, and extract name-value pairs
        default_eventToAlign = 'hitLick'; %

        p = inputParser; % create parser object
        addRequired(p, 'filePath');
        addRequired(p, 'numPeriFrames');
        addRequired(p, 'trial'); 
       
        addParameter(p, 'eventToAlign', default_eventToAlign);

        parse(p, filePath, numPeriFrames, trial, vargs{:})
    end
