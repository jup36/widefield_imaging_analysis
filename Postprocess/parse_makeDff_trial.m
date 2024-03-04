function p = parse_makeDff_trial(stack, stack_base_index, vargs)
% parse input, and extract name-value pairs
default_dffMethod = 'mean'; % default method for dff is moving average
default_type = 'dff'; % dff or fractional
default_detrend = false; % detrend or not

p = inputParser; % create parse object
addRequired(p, 'stack')
addRequired(p, 'stack_base_index')

addParameter(p, 'dffMethod', default_dffMethod)
addParameter(p, 'type', default_type)
addParameter(p, 'detrend', default_detrend)

parse(p, stack, stack_base_index, vargs{:})
end