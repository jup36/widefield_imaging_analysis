function p = parse_postprocessing(filePath, vargs)
% parse input, and extract name-value pairs
default_imageToUse = 'ome_stack'; % default image sets to use, it must be 'ome_stack' or 'dff_combined'
default_dffMethod = 'movingavg'; % default method for dff is moving average
default_bvCorrectLogic = true; % default logic for hemodynamic correction
default_tbytBaseWinF = 'all'; % default baseline window for trial-by-trial dff ('all' is to use all frames with negative timestamps relativev to the aligned event), otherwise it should be numeric, e.g., 1 for 1 s
default_movingWinF = 30; % default window for moving average dff (30 s)
default_saveTbytDff = true;
default_combineStacksOfSameSequenceLogic = true;

p = inputParser; % create parser object
addRequired(p, 'filePath')

addParameter(p, 'imageToUse', default_imageToUse)
addParameter(p, 'dffMethod', default_dffMethod)
addParameter(p, 'bvCorrectLogic', default_bvCorrectLogic)
addParameter(p, 'tbytBaseWinF', default_tbytBaseWinF)
addParameter(p, 'movingWinF', default_movingWinF)
addParameter(p, 'saveTbytDff', default_saveTbytDff)
addParameter(p, 'combineStacksOfSameSequenceLogic', default_combineStacksOfSameSequenceLogic)

parse(p, filePath, vargs{:})
end
