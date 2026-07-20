function [groupOtherBetaTable, animalOtherBetaTable, groupNames] = computeMotifOtherBetaByGroup(glmRezC, glmLabelC, groupDefs, motifId, varargin)
%COMPUTEMOTIFOTHERBETABYGROUP  Group-mean "other" (non-tone) predictor beta for one motif.
%
% SYNOPSIS
%   [groupOtherBetaTable, animalOtherBetaTable, groupNames] = ...
%       computeMotifOtherBetaByGroup(glmRezC, glmLabelC, groupDefs, motifId, ...)
%
% DESCRIPTION
%   The non-tone ("other") predictor companion to computeMotifToneBetaByGroup.m.
%   For a single, user-selected motif, collapses every NON-tone predictor
%   into one mean-beta value per "base name" (e.g. all 9 lick_rc# columns
%   -> one "lick" value; a single scalar predictor like 'pupil' -> just
%   itself) -- the exact same collapsing plot_beta_groupMeanBars.m does
%   with 'SplitOtherByBase', true, 'plotTonePredictors', false -- but now
%   averaged first across each animal's valid sessions, then across
%   animals WITHIN each group defined in groupDefs (e.g. fast/slow
%   learners), instead of being computed for one session at a time.
%
%   Base names are discovered dynamically per session (whatever
%   non-tone predictors happen to be present -- e.g. continuous
%   behavioral predictors like pupil/whisker/nosetip/locomVel are only
%   included in a session if that data was available, per the original
%   GLM design), so different animals/sessions may contribute to
%   different subsets of base names. This is handled by returning a
%   LONG-FORMAT table rather than a fixed-width one -- see OUTPUTS.
%
% INPUTS
%   glmRezC   : [nAnimals x nSessions] cell array; glmRezC{i,j} is either
%               [] or a struct with fields 'beta' ([P x K]) and 'X_names'
%               (1xP cellstr).
%   glmLabelC : cell array of label strings containing the animal ID
%               (e.g. 'm1045'), same convention as collect_cvR2_perAnimal.m.
%   groupDefs : scalar struct, field name = group label, value = cellstr
%               of animal IDs (e.g. groupDefs.fast = {...}; groupDefs.slow = {...};).
%   motifId   : scalar motif index (column of beta) to extract.
%
% NAME-VALUE ARGS
%   'excludeBases' : cellstr of base names (case-insensitive) to drop
%                    entirely (e.g. {'lick','water','airpuff'} to keep
%                    only continuous behavioral predictors). Default: {}
%                    (include every non-tone base found).
%
% OUTPUTS
%   groupOtherBetaTable : long-format table, one row per (group, baseName)
%                         pair that had at least one contributing animal.
%                         Columns: group, baseName, meanBeta, nAnimals
%                         (how many animals in that group contributed).
%   animalOtherBetaTable : long-format table, one row per (animal, baseName)
%                         pair that had at least one contributing session.
%                         Columns: animalID, group, baseName, meanBeta,
%                         nSessions (how many sessions contributed).
%                         This is the STAGE-1 (session-averaged) result --
%                         inspect this directly if a specific base/animal
%                         combination looks worth checking before trusting
%                         the group-level mean.
%   groupNames           : cellstr of group labels (fieldnames of groupDefs).
%
% EXAMPLE
%   groupDefs.fast = fast_learners;
%   groupDefs.slow = slow_learners;
%   [groupOtherBetaTable, animalOtherBetaTable, groupNames] = ...
%       computeMotifOtherBetaByGroup(glmRezC, glmLabelC, groupDefs, 20);
%
% See also: plotMotifOtherBetaByGroup, plot_beta_groupMeanBars, computeMotifToneBetaByGroup

p = inputParser;
p.addParameter('excludeBases', {}, @(x) iscell(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;
excludeBases = lower(cellstr(string(opt.excludeBases)));

[nAnimals, nSessions] = size(glmRezC);
animalIDs = extract_animal_ids_(glmLabelC, nAnimals);

% -------- Stage 1: per-animal, per-base mean across that animal's sessions --------
animalID_rows   = {};
baseName_rows   = {};
meanBeta_rows   = [];
nSessions_rows  = [];

for i = 1:nAnimals
    baseVals = containers.Map('KeyType','char','ValueType','any');   % baseName -> vector of session-level means

    for j = 1:nSessions
        s = glmRezC{i,j};
        if isempty(s) || ~isstruct(s) || ~isfield(s,'beta') || ~isfield(s,'X_names')
            continue;
        end
        if motifId > size(s.beta, 2)
            continue;
        end

        bCol = s.beta(:, motifId);
        [typeKey, baseKey] = classify_predictor_names_(s.X_names);
        mOther = (typeKey == "other");
        if ~any(mOther)
            continue;
        end

        basesThisSession = unique(baseKey(mOther), 'stable');
        for bI = 1:numel(basesThisSession)
            bName = char(basesThisSession(bI));
            if ismember(bName, excludeBases)
                continue;
            end
            mm = mOther & (baseKey == basesThisSession(bI));
            sessionMean = mean(bCol(mm), 'omitnan');

            if ~isKey(baseVals, bName)
                baseVals(bName) = sessionMean;
            else
                baseVals(bName) = [baseVals(bName), sessionMean]; %#ok<NASGU>
            end
        end
    end

    baseNamesThisAnimal = keys(baseVals);
    for bI = 1:numel(baseNamesThisAnimal)
        bName = baseNamesThisAnimal{bI};
        vals  = baseVals(bName);

        animalID_rows{end+1,1}  = animalIDs{i};      %#ok<AGROW>
        baseName_rows{end+1,1}  = bName;              %#ok<AGROW>
        meanBeta_rows(end+1,1)  = mean(vals, 'omitnan'); %#ok<AGROW>
        nSessions_rows(end+1,1) = numel(vals);        %#ok<AGROW>
    end
end

assert(~isempty(animalID_rows), ...
    'No animal had any valid session with a non-tone ("other") predictor for motif %d.', motifId);

% -------- group assignment per animal --------
groupNames = fieldnames(groupDefs);
groupTag_rows = repmat({'unassigned'}, numel(animalID_rows), 1);
for gi = 1:numel(groupNames)
    memberIDs = groupDefs.(groupNames{gi});
    groupTag_rows(ismember(animalID_rows, memberIDs)) = groupNames(gi);
end

animalOtherBetaTable = table(animalID_rows(:), groupTag_rows(:), baseName_rows(:), meanBeta_rows(:), nSessions_rows(:), ...
    'VariableNames', {'animalID','group','baseName','meanBeta','nSessions'});

% -------- Stage 2: per-group, per-base mean across that group's animals --------
uniqueBases = unique(animalOtherBetaTable.baseName);

group_rows    = {};
base_rows2    = {};
meanBeta_rows2 = [];
nAnimals_rows2 = [];

for gi = 1:numel(groupNames)
    gName = groupNames{gi};
    for bI = 1:numel(uniqueBases)
        bName = uniqueBases{bI};
        rowsI = strcmp(animalOtherBetaTable.group, gName) & strcmp(animalOtherBetaTable.baseName, bName);
        if ~any(rowsI)
            continue;
        end
        group_rows{end+1,1}     = gName; %#ok<AGROW>
        base_rows2{end+1,1}     = bName; %#ok<AGROW>
        meanBeta_rows2(end+1,1) = mean(animalOtherBetaTable.meanBeta(rowsI), 'omitnan'); %#ok<AGROW>
        nAnimals_rows2(end+1,1) = sum(rowsI); %#ok<AGROW>
    end
end

groupOtherBetaTable = table(group_rows(:), base_rows2(:), meanBeta_rows2(:), nAnimals_rows2(:), ...
    'VariableNames', {'group','baseName','meanBeta','nAnimals'});

end % function


% ===== helper: classify predictor by tone type + "other" base name =====
% (identical to plot_beta_groupMeanBars.m's local helper, duplicated for
% self-containment)
function [typeKey, baseKey] = classify_predictor_names_(nameC)
n = numel(nameC);
typeKey = repmat("other", n, 1);
baseKey = repmat("other", n, 1);

for i = 1:n
    nm = string(nameC{i});

    if contains(nm, "toneOnGo", "IgnoreCase", true)
        typeKey(i) = "tOnG";   baseKey(i) = "tOnG";
    elseif contains(nm, "toneOffGo", "IgnoreCase", true)
        typeKey(i) = "tOffG";  baseKey(i) = "tOffG";
    elseif contains(nm, "toneOnNoGo", "IgnoreCase", true)
        typeKey(i) = "tOnNg";  baseKey(i) = "tOnNg";
    elseif contains(nm, "toneOffNoGo", "IgnoreCase", true)
        typeKey(i) = "tOffNg"; baseKey(i) = "tOffNg";
    else
        base = regexprep(nm, '_rc\d+$', '');
        base = regexprep(base, '_+', '');
        baseKey(i) = lower(base);
    end
end
end


% ===== helper: pull animal ID (e.g. 'm1045') from glmLabelC row =====
function animalIDs = extract_animal_ids_(glmLabelC, nAnimals)
animalIDs = cell(nAnimals, 1);

for i = 1:nAnimals
    if size(glmLabelC,1) < i
        rowLabels = {};
    elseif size(glmLabelC,2) > 1
        rowLabels = glmLabelC(i,:);
    else
        rowLabels = glmLabelC(i,1);
    end

    idFound = '';
    for j = 1:numel(rowLabels)
        lbl = rowLabels{j};
        if (ischar(lbl) || isstring(lbl)) && ~isempty(lbl)
            tok = regexp(char(lbl), 'm\d{4}', 'match', 'once');
            if ~isempty(tok)
                idFound = tok;
                break;
            end
        end
    end

    if isempty(idFound)
        idFound = sprintf('Animal%02d', i);
    end
    animalIDs{i} = idFound;
end
end
