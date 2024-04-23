function p = parseInput_tbyt_auditory_gng_behavior( filePath, vargs )
% parse input, and extract name-value pairs
default_preToneWin = 1; % pre-Tone period for event window
default_postToneWin = 4; % post-Tone window period for event window

p = inputParser; % create parser object
addRequired(p,'filePath');
addParameter(p,'preToneWin', default_preToneWin);
addParameter(p,'postToneWin', default_postToneWin);

parse(p, filePath, vargs{:})
end