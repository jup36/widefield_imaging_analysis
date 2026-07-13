function [glmCvR2C_collect, animalIDs] = collect_cvR2_perAnimal(glmEvRezC, glmLabelC)
%COLLECT_CVR2_PERANIMAL  Collect cross-validated global R^2 across sessions, per animal.
%
% SYNOPSIS
%   [glmCvR2C_collect, animalIDs] = collect_cvR2_perAnimal(glmEvRezC, glmLabelC)
%
% DESCRIPTION
%   glmEvRezC is a [nAnimals x nSessions] cell array; each non-empty cell
%   holds a struct with field 'cvR2_global' (e.g. glmEvRezC{2,2}.cvR2_global).
%   Empty cells (missing/failed sessions) are skipped. For each animal
%   (row), the valid cvR2_global values are concatenated IN COLUMN ORDER
%   into a single vector.
%
%   Animal IDs are pulled from glmLabelC by regex-matching the pattern
%   'm####' (e.g. 'm1045') against whatever string(s) are found in that
%   animal's row -- this works whether glmLabelC is [nAnimals x nSessions]
%   (one label per session) or [nAnimals x 1] (one label per animal).
%
% INPUTS
%   glmEvRezC : [nAnimals x nSessions] cell array; glmEvRezC{i,j} is either
%               [] (missing session) or a 1x1 struct with field cvR2_global.
%   glmLabelC : cell array of label strings containing the animal ID
%               (e.g. 'm1045' or 'm1045_122424...'), same number of rows
%               as glmEvRezC (nAnimals x nSessions or nAnimals x 1).
%
% OUTPUTS
%   glmCvR2C_collect : [nAnimals x 1] cell array; glmCvR2C_collect{i} is a
%                      column vector of cvR2_global values for animal i,
%                      concatenated across its valid (non-empty) sessions
%                      in original column order. Empty cells are skipped,
%                      not NaN-padded -- so index k in this vector is the
%                      k-th VALID session for that animal, not necessarily
%                      the k-th calendar session.
%   animalIDs        : [nAnimals x 1] cellstr of animal IDs (e.g. 'm1045').
%                      Falls back to 'AnimalNN' if no match is found in
%                      glmLabelC for that row.
%
% NOTE
%   Because empty sessions are skipped rather than NaN-padded, the
%   resulting "session index" (1, 2, 3, ...) is sequential valid-session
%   count, not calendar session number. If animals have missing sessions
%   in different places, session index k does not necessarily mean the
%   same calendar timepoint across animals -- keep this in mind before
%   drawing cross-animal conclusions from x-axis position (e.g. fast vs.
%   slow learner comparisons).
%
% EXAMPLE
%   [glmCvR2C_collect, animalIDs] = collect_cvR2_perAnimal(glmEvRezC, glmLabelC);

[nAnimals, nSessions] = size(glmEvRezC);

glmCvR2C_collect = cell(nAnimals, 1);
for i = 1:nAnimals
    vals = [];
    for j = 1:nSessions
        s = glmEvRezC{i,j};
        if ~isempty(s) && isstruct(s) && isfield(s, 'cvR2_global') && ~isempty(s.cvR2_global)
            vals(end+1) = s.cvR2_global; %#ok<AGROW>
        end
    end
    glmCvR2C_collect{i} = vals(:);
end

animalIDs = extract_animal_ids_(glmLabelC, nAnimals);

end % function


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
