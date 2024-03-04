function [stack_out, avgproj] = makeDFF_trial(stack, stack_base_index, varargin)
%This function takes the imput raw stack and and makes a dff along the third dimension on a trial-by-trial basis.
p = parse_makeDff_trial(stack, stack_base_index, varargin);
%p = parse_makeDff_trial(stack, stack_base_index, {'dffMethod', 'mean', 'type', 'dff', 'detrend', false});

% choose type of averaging
if strcmp(p.Results.dffMethod, 'mean')
    avgproj  = nanmean(stack(:, :, stack_base_index),3);
elseif strcmp(p.Results.dffMethod,'median')
    avgproj  = nanmedian(stack(:, :, stack_base_index),3);
else
    error('unknown dff method!');
end

% calculate df and dff, note that 'dff_simple =
% (stack-avgproj)./avgproj*100;' works as well in recent MATLAB versions (2016 or later)
% which support implicit expansion. 
df = bsxfun(@minus, stack, avgproj);
%dff_simple = (stack-avgproj)./avgproj*100;

if strcmp(p.Results.type, 'dff')
   stack_out = bsxfun(@rdivide, df, avgproj).*100;
elseif strcmp(p.Results.type, 'fractional')
   stack_out = bsxfun(@rdivide, df, avgproj); 
else
    error('unknown type of calculation (use dff or fractional)!')
end

%linear detrend
if p.Results.detrend
    stack_out = detrendNaN3(stack_out);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

end