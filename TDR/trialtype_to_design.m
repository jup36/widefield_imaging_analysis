function [X, names, info] = trialtype_to_design(trialType, varargin)
%TRIALTYPE_TO_DESIGN  Build per-trial dummy regressors for trial type.
%
% [X, names, info] = trialtype_to_design(trialType, ... )
%
% INPUT
%   trialType : N×1 categorical or cellstr/char array with labels
%               (default universe: {'Hit','Miss','CR','FA'})
%
% NAME–VALUE OPTIONS
%   'catsAll'   : cellstr of all categories in desired order
%                 (default {'Hit','Miss','CR','FA'})
%   'ref'       : reference category to drop (default 'Hit')
%   'dropEmpty' : logical; drop columns for categories absent this session (default true)
%   'prefix'    : char/str; column name prefix (default 'TT_')
%   'warn'      : logical; print info about dropped/absent (default true)
%
% OUTPUT
%   X      : N×P dummy matrix (P = #non-reference categories present)
%   names  : 1×P cellstr of column names
%   info   : struct with fields:
%              .catsAll, .ref, .keptCats, .droppedCats, .absentCats
%
% Notes
%   - Reference coding: each column is an indicator for a non-reference category.
%   - If only the reference category is present, X can be empty.
%
% Junchol Park / Buschman Lab — 2025

% ---------- parse args ----------
p = inputParser;
p.addParameter('catsAll',   {'Hit','Miss','CR','FA'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('ref',       'Hit',                    @(s) ischar(s) || isstring(s));
p.addParameter('dropEmpty', true,                     @(v) islogical(v) || ismember(v,[0 1]));
p.addParameter('prefix',    'TT_',                    @(s) ischar(s) || isstring(s));
p.addParameter('warn',      true,                     @(v) islogical(v) || ismember(v,[0 1]));
p.parse(varargin{:});
opt = p.Results;

catsAll = cellstr(string(opt.catsAll));
ref     = char(string(opt.ref));
if ~ismember(ref, catsAll)
    error('Reference "%s" is not in catsAll.', ref);
end

% ---------- coerce to categorical with fixed order ----------
yCat = categorical(trialType, catsAll, catsAll);
N = numel(yCat);

% ---------- build reference-coded dummies ----------
nonRefCats = catsAll(~strcmp(catsAll, ref));
D = zeros(N, numel(nonRefCats));
for k = 1:numel(nonRefCats)
    D(:,k) = double(yCat == nonRefCats{k});
end
names = strcat(opt.prefix, nonRefCats);

% ---------- drop empty categories if requested ----------
isEmpty = ~any(D~=0, 1);
absentCats = nonRefCats(isEmpty);
if opt.dropEmpty && any(isEmpty)
    D(:, isEmpty) = [];
    names(isEmpty) = [];
end
keptCats    = nonRefCats(~isEmpty);
droppedCats = ternary(opt.dropEmpty, absentCats, {});

X = D;

% ---------- info & optional logging ----------
info = struct( ...
    'catsAll',     {catsAll}, ...
    'ref',         ref, ...
    'keptCats',    {keptCats}, ...
    'droppedCats', {droppedCats}, ...
    'absentCats',  {absentCats});

if opt.warn
    if isempty(keptCats)
        fprintf('[trialtype_to_design] Only reference category "%s" present; X is empty.\n', ref);
    end
    if ~isempty(absentCats)
        if opt.dropEmpty
            fprintf('[trialtype_to_design] Dropped absent cats: %s\n', strjoin(absentCats, ', '));
        else
            fprintf('[trialtype_to_design] Absent cats (kept as zero cols): %s\n', strjoin(absentCats, ', '));
        end
    end
end
end

% --- tiny local helper ---
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
