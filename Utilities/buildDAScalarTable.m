function T = buildDAScalarTable(rezC, sessT, varargin)
%BUILDDASCALARTABLE
%   Long table of session-level mean DA, one row per session x trial stream
%   x time window, with d', learner group, day-4 status and animal. Shared
%   by the d'-correlation plots so the scalar is defined in exactly one
%   place -- two copies of this logic would drift.
%
%   The DA value is the session's mean PETH averaged over the window, which
%   equals the mean over trials of each trial's window mean (both linear).
%
%   T = buildDAScalarTable(rezC, sessT, ...)
%
% NAME-VALUE
%   'signal'       : 'global' (default) or motif number 1..K
%   'trialTypes'   : {'hit','cr'}
%   'windows'      : [0 2; 2 4]
%   'windowNames'  : {'cue','response'}
%   'fastLearners' : {'m1044','m1045','m1092','m1094'}
%   'slowLearners' : {'m1048','m1049','m1613','m1859','m1873'}
%
% Animals in neither group (e.g. m1237) are omitted. Sessions with NaN d'
% are KEPT here and filtered by the caller, so the table stays a complete
% record of what was computed.

p = inputParser;
p.addParameter('signal', 'global', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && isscalar(x)));
p.addParameter('trialTypes', {'hit','cr'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('windows', [0 2; 2 4], @(x) isnumeric(x) && size(x,2) == 2);
p.addParameter('windowNames', {'cue','response'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('fastLearners', {'m1044','m1045','m1092','m1094'}, @iscell);
p.addParameter('slowLearners', {'m1048','m1049','m1613','m1859','m1873'}, @iscell);
p.parse(varargin{:});
opt = p.Results;

trialTypes  = cellstr(opt.trialTypes);
windowNames = cellstr(opt.windowNames);
W = opt.windows;
assert(size(W,1) == numel(windowNames), 'windows and windowNames must match in number.');
assert(ismember('isDay4', sessT.Properties.VariableNames), ...
    'sessT has no isDay4 column -- rerun batch_DAPeth_globalAndMotifWeighted with day4MarkC.');

if isnumeric(opt.signal)
    k = opt.signal;
    getTrace = @(pe) pe.motifMean(k, :);
else
    assert(strcmpi(opt.signal, 'global'), 'signal must be ''global'' or a motif number.');
    getTrace = @(pe) pe.globalMean;
end

tint = rezC{sessT.a(1), sessT.s(1)}.meta.tint;
winI = arrayfun(@(i) tint >= W(i,1) & tint < W(i,2), 1:size(W,1), 'UniformOutput', false);
for i = 1:numel(winI)
    assert(any(winI{i}), 'Window %s selects no bins of the time grid.', windowNames{i});
    % A window reaching past the grid is silently truncated to what exists,
    % which would make it shorter than its label says -- flag it.
    if W(i,1) < tint(1) || W(i,2) > tint(end) + (tint(2) - tint(1))/2
        warning('DAScalar:WindowPastGrid', ...
            ['Window %s [%g %g] extends past the time grid [%g %g]; only the ' ...
             'overlapping %g-%g s is averaged.'], windowNames{i}, W(i,1), W(i,2), ...
            tint(1), tint(end), max(W(i,1), tint(1)), min(W(i,2), tint(end)));
    end
end

animal = {}; group = {}; hdr = {}; tt_ = {}; win_ = {}; dp = []; da = []; d4 = []; so = [];
for i = 1:height(sessT)
    r = rezC{sessT.a(i), sessT.s(i)};
    aId = char(sessT.animal{i});
    if ismember(aId, opt.fastLearners), g = 'fast';
    elseif ismember(aId, opt.slowLearners), g = 'slow';
    else, continue;
    end
    for t = 1:numel(trialTypes)
        tt = trialTypes{t};
        if ~isfield(r.peth, tt) || r.peth.(tt).n == 0, continue; end
        y = getTrace(r.peth.(tt));
        for w = 1:size(W,1)
            animal{end+1,1} = aId;                          %#ok<AGROW>
            group{end+1,1}  = g;                            %#ok<AGROW>
            hdr{end+1,1}    = char(sessT.header{i});        %#ok<AGROW>
            tt_{end+1,1}    = tt;                           %#ok<AGROW>
            win_{end+1,1}   = windowNames{w};               %#ok<AGROW>
            dp(end+1,1)     = sessT.dprime(i);              %#ok<AGROW>
            da(end+1,1)     = mean(y(winI{w}), 'omitnan');  %#ok<AGROW>
            d4(end+1,1)     = sessT.isDay4(i);              %#ok<AGROW>
            so(end+1,1)     = sessT.sessOrder(i);           %#ok<AGROW>
        end
    end
end

T = table(animal, group, hdr, tt_, win_, dp, da, logical(d4), so, ...
    'VariableNames', {'animal','group','header','trialType','window', ...
                      'dprime','DA','isDay4','sessOrder'});
end