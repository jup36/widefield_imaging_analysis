%% checkScottyDependencies.m
%
% Checks whether each file required by
% motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func (as
% reported by matlab.codetools.requiredFilesAndProducts) already exists
% somewhere on Scotty, searched by filename only (not exact path, since
% local Windows folder structure does not need to match the bucket
% structure -- MATLAB path resolution is flat).
%
% Requires an active s_conn (see earlier ssh2_command_scotty('connect', ...)).
%
% USAGE
%   [fList, ~] = matlab.codetools.requiredFilesAndProducts( ...
%       'motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func.m');
%   checkScottyDependencies(fList, s_conn);

function checkScottyDependencies(fList, s_conn)

% Scoped search root -- fast, since it's this lab's own subtree rather
% than a filesystem-wide crawl (which is what hung earlier when searching
% from /).
searchRoot = '/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis';

fprintf('\n============================================================\n');
fprintf('Checking %d dependencies against Scotty (%s)\n', numel(fList), searchRoot);
fprintf('============================================================\n');

missing = {};

for i = 1:numel(fList)
    [~, fname, ext] = fileparts(fList{i});
    targetName = [fname ext];

    checkCmd = sprintf('find -L "%s" -iname "%s" 2>/dev/null', searchRoot, targetName);
    resp = ssh2_command_scotty(s_conn, checkCmd);

    resultLines = resp.command_result;
    if ~iscell(resultLines)
        resultLines = cellstr(resultLines);
    end
    % Drop the SSH banner/MOTD line the same way the submitter does --
    % anything not looking like an actual filesystem path.
    resultLines = resultLines(startsWith(strtrim(resultLines), '/'));

    if isempty(resultLines)
        status = 'MISSING';
        missing{end+1} = fList{i}; %#ok<AGROW>
    else
        status = 'FOUND';
    end

    fprintf('%-32s %-8s', targetName, status);
    if ~isempty(resultLines)
        fprintf('  -> %s', strjoin(resultLines, ', '));
    end
    fprintf('\n');
end

fprintf('\n============================================================\n');
if isempty(missing)
    fprintf('All dependencies found on Scotty.\n');
else
    fprintf('%d dependency file(s) MISSING on Scotty -- copy these from Windows:\n', numel(missing));
    for i = 1:numel(missing)
        fprintf('  %s\n', missing{i});
    end
end
fprintf('============================================================\n');

end