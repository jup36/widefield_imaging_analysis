function cmapOut = animalColor(selectMice, mIdC, cmapRef)
% animalColor
%   Return RGB colors for a selected list of animals based on a reference
%   animal list and its colormap.
%
%   cmapOut = animalColor(selectMice, mIdC, cmapRef)
%
% INPUTS
%   selectMice : cell array or string array of animal IDs to plot
%                e.g. {'m1044','m1045','m1092'}
%
%   mIdC       : reference list of animal IDs (cell or string)
%                e.g. {'m1237','m1092','m1094',...}
%
%   cmapRef    : [N x 3] RGB colormap corresponding to mIdC
%                e.g. cmap10
%
% OUTPUT
%   cmapOut    : [numel(selectMice) x 3] RGB colormap
%                rows correspond to selectMice order
%
% EXAMPLE
%   cmapSel = animalColor({'m1044','m1092'}, mIdC, cmap10);
%
% NOTES
%   • Matching is exact string matching
%   • Errors if a requested animal is not found in mIdC
%

% -------------------- input normalization --------------------
selectMice = string(selectMice(:));
mIdC       = string(mIdC(:));

if size(cmapRef,1) ~= numel(mIdC)
    error('animalColor:SizeMismatch', ...
        'cmapRef must have one row per entry in mIdC.');
end

% -------------------- match and extract colors --------------------
nSel = numel(selectMice);
cmapOut = nan(nSel, 3);

for i = 1:nSel
    idx = find(mIdC == selectMice(i), 1, 'first');

    if isempty(idx)
        error('animalColor:AnimalNotFound', ...
            'Animal "%s" not found in reference mIdC.', selectMice(i));
    end

    cmapOut(i, :) = cmapRef(idx, :);
end

end
