function [headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_(headerC, glmRezC, trIdC)
% Trim trailing columns that are entirely empty ([]) in whichever arrays have them.
% This prevents size-mismatch errors when one input has an extra empty column at the end.

assert(iscell(headerC) && iscell(glmRezC) && iscell(trIdC), 'Inputs must be cell arrays.');
[J1,S1] = size(headerC);
[J2,S2] = size(glmRezC);
[J3,S3] = size(trIdC);

% If row counts disagree, don't guess.
if ~(J1==J2 && J1==J3)
    error('Row counts differ: headerC=%dx%d, glmRezC=%dx%d, trIdC=%dx%d', J1,S1,J2,S2,J3,S3);
end

Smax = max([S1 S2 S3]);

% Pad smaller ones with empty so we can test columns uniformly
if S1 < Smax, headerC(:,end+1:Smax) = {[]}; end
if S2 < Smax, glmRezC(:,end+1:Smax) = {[]}; end
if S3 < Smax, trIdC(:,end+1:Smax)   = {[]}; end

% Find last column that has ANY non-empty content in ANY of the three inputs
keepLast = 0;
for s = 1:Smax
    anyNonEmpty = any(~cellfun(@isempty, headerC(:,s))) || ...
                  any(~cellfun(@isempty, glmRezC(:,s)))  || ...
                  any(~cellfun(@isempty, trIdC(:,s)));
    if anyNonEmpty
        keepLast = s;
    end
end

if keepLast == 0
    % all empty everywhere; keep a 0-column shape consistent with original rows
    headerC = cell(J1,0);
    glmRezC = cell(J1,0);
    trIdC   = cell(J1,0);
else
    headerC = headerC(:,1:keepLast);
    glmRezC = glmRezC(:,1:keepLast);
    trIdC   = trIdC(:,1:keepLast);
end
end