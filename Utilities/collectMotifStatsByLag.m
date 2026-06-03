function [stats_trainPevC, stats_trainNmotifC, ...
          stats_testPevC, stats_testNmotifC] = ...
          collectMotifStatsByLag(filePath_base, mlist, lagNum, chunkSearchStr)

% collectMotifStatsByLag
%
% Collects motif-fitting statistics across mice/sessions/chunks for a given
% lag condition and chunk-file search string.
%
% Inputs:
%   filePath_base
%       Base directory containing preprocessed motif folders.
%
%   mlist
%       Cell array of mouse IDs, e.g. {'m1044', 'm1045'}.
%
%   lagNum
%       Motif lag condition.
%       Use 1, 3, 5 for explicit lag folders.
%       Use 10 for original/default motif folders ending in 'motif' and
%       not containing 'lag'.
%
%   chunkSearchStr
%       Search pattern for chunk files inside each session folder.
%       Examples:
%           '*green*chunk*.mat'
%           '*red*chunk*.mat'
%
% Outputs:
%   stats_trainPevC
%   stats_trainNmotifC
%   stats_testPevC
%   stats_testNmotifC
%
% Each output is a mouse × session cell array. Each cell contains the mean
% value across chunks for that mouse/session, or NaN if no valid chunks exist.

stats_trainPevC    = cell(numel(mlist), 20); 
stats_trainNmotifC = cell(numel(mlist), 20); 
stats_testPevC     = cell(numel(mlist), 20);  
stats_testNmotifC  = cell(numel(mlist), 20); 

for f = 1:numel(mlist)

%     %% Find motif session folders
%     if lagNum == 10
% 
%         % Original/default motifs:
%         % folder name should end with 'motif' and should not contain 'lag'
%         sessions = GrabFiles_sort_trials([mlist{f} '*motif'], 0, {filePath_base});  
% 
%         [~, sessionNames] = cellfun(@fileparts, sessions, 'UniformOutput', false);
% 
%         keepI = endsWith(sessionNames, 'motif') & ...
%                 ~contains(sessionNames, 'lag', 'IgnoreCase', true);
% 
%         sessions = sessions(keepI);
% 
%     else

        % Explicit lag folders
        sessions = GrabFiles_sort_trials( ...
            [mlist{f} '*motif*lag' num2str(lagNum)], ...
            0, {filePath_base});  

%     end

    fprintf('\nMouse %s | Lag %d | Search: %s | %d sessions found\n', ...
        mlist{f}, lagNum, chunkSearchStr, numel(sessions));

    %% Iterate through sessions
    for ff = 1:numel(sessions)

        % Initialize as NaN by default
        stats_trainPevC{f, ff}    = NaN;
        stats_trainNmotifC{f, ff} = NaN;
        stats_testPevC{f, ff}     = NaN;
        stats_testNmotifC{f, ff}  = NaN;

        header = extract_date_animalID_header(sessions{ff}); 

        % Combine session header with user-specified chunk search string.
        %
        % Example:
        %   header = 'm1045_122424'
        %   chunkSearchStr = '*green*chunk*.mat'
        %   searchStr = 'm1045_122424*green*chunk*.mat'
        searchStr = [header, chunkSearchStr]; 

        filePath_chunk = GrabFiles_sort_trials(searchStr, 0, sessions(ff));

        if isempty(filePath_chunk)
            fprintf('%s | lag %d | session %d/%d: no chunks found with %s\n', ...
                mlist{f}, lagNum, ff, numel(sessions), searchStr);
            continue
        end

        %% Initialize temporary numeric vectors
        train_pev    = nan(numel(filePath_chunk), 1);
        train_nmotif = nan(numel(filePath_chunk), 1);
        test_pev     = nan(numel(filePath_chunk), 1);
        test_nmotif  = nan(numel(filePath_chunk), 1);

        %% Load each chunk
        for j = 1:numel(filePath_chunk)

            try
                temp = load(filePath_chunk{j}, 'stats_train', 'stats_test');
            catch ME
                warning('Could not load chunk file:\n%s\nError: %s', ...
                    filePath_chunk{j}, ME.message);
                continue
            end

            if isfield(temp, 'stats_train') && isfield(temp.stats_train, 'pev')
                train_pev(j) = temp.stats_train.pev;
            end

            if isfield(temp, 'stats_train') && isfield(temp.stats_train, 'n_motifs')
                train_nmotif(j) = temp.stats_train.n_motifs;
            end

            if isfield(temp, 'stats_test') && isfield(temp.stats_test, 'pev')
                test_pev(j) = temp.stats_test.pev;
            end

            if isfield(temp, 'stats_test') && isfield(temp.stats_test, 'n_motifs')
                test_nmotif(j) = temp.stats_test.n_motifs;
            end

        end

        %% Average across chunks within this session
        stats_trainPevC{f, ff}    = mean(train_pev, 'omitnan');
        stats_trainNmotifC{f, ff} = mean(train_nmotif, 'omitnan');

        stats_testPevC{f, ff}     = mean(test_pev, 'omitnan');
        stats_testNmotifC{f, ff}  = mean(test_nmotif, 'omitnan');

        fprintf('%s | lag %d | session %d/%d | %d chunks | train PEV %.3f | test PEV %.3f | train nMotif %.2f | test nMotif %.2f\n', ...
            mlist{f}, lagNum, ff, numel(sessions), numel(filePath_chunk), ...
            stats_trainPevC{f, ff}, stats_testPevC{f, ff}, ...
            stats_trainNmotifC{f, ff}, stats_testNmotifC{f, ff});

    end
end
end
