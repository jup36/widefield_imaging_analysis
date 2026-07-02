%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sessionID = localMakeBaselineSessionIDFromLKFolder(folderPath)
% localMakeBaselineSessionIDFromLKFolder
%
% Converts a baseline LK folder basename into a session ID.
%
% Example:
%   Input folder:
%       m1048_122424_base_day4-5_img_motif_lag1_k10
%
%   Output sessionID:
%       m1048_122424_base_day4-5

[~, folderName] = fileparts(folderPath);

sessionID = regexprep(folderName, '_img_motif_lag\d+_k\d+.*$', '');

if isempty(sessionID)
    sessionMatch = regexp(folderName, 'm\d{4}_\d{6}_base[^_]*.*?(?=_img_motif)', ...
        'match', 'once');
    
    if ~isempty(sessionMatch)
        sessionID = sessionMatch;
    else
        sessionID = folderName;
    end
end

end