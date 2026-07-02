%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function firstPath = localFirstCell(pathC)
% localFirstCell
%
% Safely returns the first element of a cell array of paths.
% Returns '' if input is empty.
%
% This avoids fragile cell2mat(...) behavior when find_keyword_containing_folder
% returns zero or multiple matches.

if isempty(pathC)
    firstPath = '';
elseif iscell(pathC)
    firstPath = pathC{1};
else
    firstPath = pathC;
end

end