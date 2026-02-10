function glmTargetCols = glmNameToColumns(glmNameListC, targetNameC)
%GLMNAMETOCOLUMNS  Map target GLM axis names to column indices.
%
% glmTargetCols = glmNameToColumns(glmNameListC, targetNameC)
%
% INPUTS
%   glmNameListC : 1xK (or Kx1) cell array of GLM axis names
%                  e.g. {'GoToneOn_1','GoToneOn_2',...}
%   targetNameC  : cell array (or string array) of names to find
%                  e.g. {'GoToneOn_1','GoToneOn_2','GoToneOn_3'}
%
% OUTPUT
%   glmTargetCols : numeric row vector of column indices into glmNameListC
%                   (exact-match only; order follows targetNameC)
%
% NOTES
%   - Exact string match only (case-sensitive)
%   - Missing target names are silently ignored
%   - Safe with char/string mixtures

% -------------------- sanity --------------------
if isempty(glmNameListC) || isempty(targetNameC)
    glmTargetCols = [];
    return;
end

% normalize to string row vectors
glmNames = string(glmNameListC(:))';    % 1 x K
targets  = string(targetNameC(:))';     % 1 x T

glmTargetCols = [];

% -------------------- exact matching --------------------
for i = 1:numel(targets)
    idx = find(glmNames == targets(i), 1, 'first');
    if ~isempty(idx)
        glmTargetCols(end+1) = idx; %#ok<AGROW>
    end
end

end
