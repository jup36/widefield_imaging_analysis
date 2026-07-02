function sessionID = localMakeSessionIDFromFolder(folderPath, sourceType)
% localMakeSessionIDFromFolder
%
% For LKcombo:
%   folder basename is usually m####_######.
%
% For legacy:
%   folder basename is usually m####_######_task_dayX_img_motif.
%   We remove _img_motif and keep the full session/day label:
%       m1613_050725_task_day4-7
%
% This avoids duplicate-date collapse.

[~, folderName] = fileparts(folderPath);

sessionID = folderName;

if strcmpi(sourceType, 'legacy')

    sessionID = regexprep(sessionID, '_img_motif.*$', '');

elseif strcmpi(sourceType, 'LKcombo')

    % Usually already clean. Keep as basename.
    sessionID = folderName;

else

    % Fallback: keep basename.
    sessionID = folderName;
end

% Final fallback if something strange happens
if isempty(sessionID)
    sessionMatch = regexp(folderPath, 'm\d{4}_\d{6}', 'match', 'once');

    if ~isempty(sessionMatch)
        sessionID = sessionMatch;
    else
        sessionID = folderName;
    end
end

end