function [synthBetaByGroup, synthNames, animalToneTable, groupNames, nAnimalsByGroup] = computeMotifToneBetaByGroup(glmRezC, glmLabelC, groupDefs, motifId, varargin)
%COMPUTEMOTIFTONEBETABYGROUP  Group-mean tone-predictor beta profile for one motif.
%
% SYNOPSIS
%   [synthBetaByGroup, synthNames, animalToneTable, groupNames, nAnimalsByGroup] = ...
%       computeMotifToneBetaByGroup(glmRezC, glmLabelC, groupDefs, motifId, ...)
%
% DESCRIPTION
%   For a single, user-selected motif, extracts the 4 tone-predictor
%   kernels (toneOnGo, toneOnNoGo, toneOffGo, toneOffNoGo -- each a
%   length-nBins vector of basis-coefficient betas) from EVERY session,
%   averages them across each animal's valid sessions, then averages
%   those animal-level profiles across animals WITHIN each group defined
%   in groupDefs (e.g. fast/slow learners).
%
%   The group-mean profiles are packaged as a synthetic beta column +
%   matching X_names list per group, in exactly the format
%   plot_beta_timeBinPaired.m expects -- so that existing plotting
%   function can be reused as-is (via motifId=1 on the synthetic data)
%   rather than reimplementing its layout/coloring logic.
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
%   'nBinsExpected' : expected number of raised-cosine basis functions per
%                     tone-predictor type. Sessions with a different count
%                     are skipped (with a warning) rather than silently
%                     misaligned. Default: 9.
%
% OUTPUTS
%   synthBetaByGroup : scalar struct, field name = group label, value =
%                      a [4*nBinsExpected x 1] synthetic beta vector,
%                      ordered [toneOnGo; toneOnNoGo; toneOffGo; toneOffNoGo]
%                      (each a block of nBinsExpected), ready to feed into
%                      plot_beta_timeBinPaired.m as a single-column "beta"
%                      matrix. Empty if that group had no valid animals.
%   synthNames       : 1x(4*nBinsExpected) cellstr of matching predictor
%                      names (e.g. 'toneOnGo_rc1', ...), shared across all
%                      groups.
%   animalToneTable  : table, one row per animal, columns: animalID,
%                      group, tOnG/tOnNg/tOffG/tOffNg (each an
%                      nAnimals x nBinsExpected matrix column; rows with
%                      no valid session are all-NaN).
%   groupNames       : cellstr of group labels (fieldnames of groupDefs).
%   nAnimalsByGroup  : scalar struct, field name = group label, value =
%                      number of animals in that group with valid data
%                      for this motif (for titling/n-reporting).
%
% NOTE
%   A session only contributes if it has usable values for ALL FOUR tone-
%   predictor types with a consistent, matching basis-index set (same
%   check plot_beta_timeBinPaired.m applies) -- sessions failing this are
%   skipped for that animal rather than partially included.
%
% EXAMPLE
%   groupDefs.fast = fast_learners;
%   groupDefs.slow = slow_learners;
%   [synthBetaByGroup, synthNames, animalToneTable, groupNames, nAnimalsByGroup] = ...
%       computeMotifToneBetaByGroup(glmRezC, glmLabelC, groupDefs, 20);
%
% See also: plotMotifToneBetaByGroup, plot_beta_timeBinPaired

p = inputParser;
p.addParameter('nBinsExpected', 9, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(varargin{:});
opt = p.Results;
nBins = opt.nBinsExpected;

[nAnimals, nSessions] = size(glmRezC);
animalIDs = extract_animal_ids_(glmLabelC, nAnimals);

tOnG_animal   = nan(nAnimals, nBins);
tOnNg_animal  = nan(nAnimals, nBins);
tOffG_animal  = nan(nAnimals, nBins);
tOffNg_animal = nan(nAnimals, nBins);

for i = 1:nAnimals
    stackOnG = []; stackOnNg = []; stackOffG = []; stackOffNg = [];

    for j = 1:nSessions
        s = glmRezC{i,j};
        if isempty(s) || ~isstruct(s) || ~isfield(s,'beta') || ~isfield(s,'X_names')
            continue;
        end
        if motifId > size(s.beta, 2)
            continue;
        end

        bCol = s.beta(:, motifId);
        [typeKey, idxNum] = classify_tone_predictor_(s.X_names);

        [vOnG,  okOnG]  = extract_sorted_(bCol, typeKey, idxNum, "tOnG");
        [vOnNg, okOnNg] = extract_sorted_(bCol, typeKey, idxNum, "tOnNg");
        [vOffG, okOffG] = extract_sorted_(bCol, typeKey, idxNum, "tOffG");
        [vOffNg,okOffNg]= extract_sorted_(bCol, typeKey, idxNum, "tOffNg");

        if ~(okOnG && okOnNg && okOffG && okOffNg)
            continue;   % missing/duplicate-index tone predictor type this session -- skip
        end
        if numel(vOnG)~=nBins || numel(vOnNg)~=nBins || numel(vOffG)~=nBins || numel(vOffNg)~=nBins
            warning('computeMotifToneBetaByGroup:binCountMismatch', ...
                'Animal %s, session %d: bin count does not match nBinsExpected=%d -- skipping this session.', ...
                animalIDs{i}, j, nBins);
            continue;
        end

        stackOnG  = [stackOnG;  vOnG(:)'];   %#ok<AGROW>
        stackOnNg = [stackOnNg; vOnNg(:)'];  %#ok<AGROW>
        stackOffG = [stackOffG; vOffG(:)'];  %#ok<AGROW>
        stackOffNg= [stackOffNg;vOffNg(:)']; %#ok<AGROW>
    end

    if ~isempty(stackOnG)
        tOnG_animal(i,:)   = mean(stackOnG,  1);
        tOnNg_animal(i,:)  = mean(stackOnNg, 1);
        tOffG_animal(i,:)  = mean(stackOffG, 1);
        tOffNg_animal(i,:) = mean(stackOffNg,1);
    end
end

groupNames = fieldnames(groupDefs);
groupCol = repmat({'unassigned'}, nAnimals, 1);
for gi = 1:numel(groupNames)
    memberIDs = groupDefs.(groupNames{gi});
    groupCol(ismember(animalIDs, memberIDs)) = groupNames(gi);
end

animalToneTable = table(animalIDs(:), groupCol(:), tOnG_animal, tOnNg_animal, tOffG_animal, tOffNg_animal, ...
    'VariableNames', {'animalID','group','tOnG','tOnNg','tOffG','tOffNg'});

% -------- group-level means + synthetic beta/name vectors --------
synthBetaByGroup = struct();
nAnimalsByGroup  = struct();

for gi = 1:numel(groupNames)
    gName = groupNames{gi};
    rowsI = strcmp(animalToneTable.group, gName) & ~isnan(animalToneTable.tOnG(:,1));
    nAnimalsByGroup.(gName) = sum(rowsI);

    if ~any(rowsI)
        warning('computeMotifToneBetaByGroup:emptyGroup', ...
            'Group "%s" has no animals with valid data for motif %d.', gName, motifId);
        synthBetaByGroup.(gName) = [];
        continue;
    end

    mOnG   = mean(animalToneTable.tOnG(rowsI,:),   1);
    mOnNg  = mean(animalToneTable.tOnNg(rowsI,:),  1);
    mOffG  = mean(animalToneTable.tOffG(rowsI,:),  1);
    mOffNg = mean(animalToneTable.tOffNg(rowsI,:), 1);

    synthBetaByGroup.(gName) = [mOnG(:); mOnNg(:); mOffG(:); mOffNg(:)];
end

synthNames = [ ...
    arrayfun(@(b) sprintf('toneOnGo_rc%d', b),    1:nBins, 'uni', 0), ...
    arrayfun(@(b) sprintf('toneOnNoGo_rc%d', b),  1:nBins, 'uni', 0), ...
    arrayfun(@(b) sprintf('toneOffGo_rc%d', b),   1:nBins, 'uni', 0), ...
    arrayfun(@(b) sprintf('toneOffNoGo_rc%d', b), 1:nBins, 'uni', 0) ...
    ];

end % function


% ===== helper: pull a type's values out sorted by basis index; flags dup/missing =====
function [v, ok] = extract_sorted_(bCol, typeKey, idxNum, targetType)
mask = (typeKey == targetType) & ~isnan(idxNum);
idx = idxNum(mask);
val = bCol(mask);
[idxSorted, ord] = sort(idx);
v = val(ord);
ok = ~isempty(v) && numel(unique(idxSorted))==numel(idxSorted);
end


% ===== helper: classify predictor by tone type + basis (time-bin) index =====
% (same logic as in plot_beta_timeBinPaired.m, duplicated for self-containment)
function [typeKey, idxNum] = classify_tone_predictor_(nameC)
n = numel(nameC);
typeKey = repmat("other", n, 1);
idxNum  = nan(n,1);

for i = 1:n
    nm = string(nameC{i});

    if contains(nm, "toneOnGo", "IgnoreCase", true)
        typeKey(i) = "tOnG";
    elseif contains(nm, "toneOffGo", "IgnoreCase", true)
        typeKey(i) = "tOffG";
    elseif contains(nm, "toneOnNoGo", "IgnoreCase", true)
        typeKey(i) = "tOnNg";
    elseif contains(nm, "toneOffNoGo", "IgnoreCase", true)
        typeKey(i) = "tOffNg";
    end

    tok = regexp(nm, '_rc(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        idxNum(i) = str2double(tok{1});
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
