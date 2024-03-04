function p = parse_stackProcessing(filePath, stack, tbytEvt, options, cmosExp, vargs)
% parse input, and extract name-value pairs
default_dffMethod = 'movingavg';    % default method for dff is moving average
default_bvCorrectLogic = true; % default logic for hemodynamic correction
default_tbytBaseWinF = 'all'; % default baseline window for trial-by-trial dff ('all' is to use all frames with negative timestamps relativev to the aligned event), otherwise it should be numeric, e.g., 1 for 1 s
default_movingWinF = 30; % default window for moving average dff (30 s)
default_saveTbytDff = true;

p = inputParser; % create parser object
addRequired(p, 'filePath')
addRequired(p, 'stack')
addRequired(p, 'options')
addRequired(p, 'tbytEvt')
addRequired(p, 'cmosExp')

addParameter(p, 'dffMethod', default_dffMethod)
addParameter(p, 'bvCorrectLogic', default_bvCorrectLogic)
addParameter(p, 'tbytBaseWinF', default_tbytBaseWinF)
addParameter(p, 'movingWinF', default_movingWinF)
addParameter(p, 'saveTbytDff', default_saveTbytDff)

parse(p, filePath, stack, tbytEvt, options, cmosExp, vargs{:})
end