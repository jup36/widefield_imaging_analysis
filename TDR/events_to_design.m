function X = events_to_design(tbytDat, winBounds, eventField, varargin)
%EVENTS_TO_DESIGN  Convert per-trial event times into a binned design matrix.
%
% X = events_to_design(tbytDat, winBounds, eventField, ...)
%
% INPUTS
%   tbytDat    : 1xN struct per trial (eventField may be [] for some trials)
%   winBounds  : nW x 2 matrix of [lb ub] edges (seconds), half-open [lb, ub)
%   eventField : char/str, name of field in tbytDat(n) containing event times
%
% NAME–VALUE PAIRS
%   'absTime'    : logical (default false). If true, treat event times as absolute
%                  and subtract tbytDat(n).evtOn to get trial-relative.
%   'trialLogic' : logical scalar or N×1 (default: true), mask selecting trials
%                  to include; masked trials remain all-zeros.
%   'mode'       : 'binary' (any event => 1) or 'count' (# events per bin). Default 'binary'.
%
% OUTPUT
%   X : N x nW matrix (binary or counts), ready to be convolved or used directly.
%
% Junchol Park / Buschman Lab — 2025

% ---------- basics ----------
N  = numel(tbytDat);
if size(winBounds,2) ~= 2
    error('winBounds must be nW x 2 [lb ub].');
end
edges = [winBounds(:,1).', winBounds(end,2)];   % nW+1 vector of edges
edges = edges(:).';                              % row
nW = size(winBounds,1);

% ---------- parse ----------
p = inputParser;
p.addParameter('absTime',    false, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('trialLogic', true(N,1), @(v)islogical(v) || isnumeric(v));
p.addParameter('mode',       'binary', @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
prm = p.Results;

% normalize trialLogic: allow scalar and numeric 0/1
mask = logical(prm.trialLogic);
if isscalar(mask), mask = repmat(mask, N, 1); end
if numel(mask) ~= N, error('trialLogic must be scalar or N×1.'); end

binaryMode = startsWith(lower(string(prm.mode)),'bin');

% ---------- allocate ----------
X = zeros(N, nW);

% ---------- loop trials ----------
for n = 1:N

    if ~mask(n)         % excluded trial -> leave zeros
        continue
    end

    if ~isfield(tbytDat, eventField) || isempty(tbytDat(n).(eventField))
        continue        % no events for this trial
    end

    evtTimes = double(tbytDat(n).(eventField));
    evtTimes = evtTimes(isfinite(evtTimes));  % drop NaNs/Inf
    if isempty(evtTimes), continue; end

    if prm.absTime
        if isfield(tbytDat, 'evtOn') && isfinite(tbytDat(n).evtOn)
            evtTimes = evtTimes - double(tbytDat(n).evtOn);
        else
            warning('events_to_design:MissingEvtOn', ...
               'Trial %d missing evtOn while absTime==true; leaving events unshifted.', n);
        end
    end

    % Bin counts using shared edges
    counts = histcounts(evtTimes, edges);   % 1 x nW

    if binaryMode
        X(n,:) = counts > 0;
    else
        X(n,:) = counts;
    end
end
end
