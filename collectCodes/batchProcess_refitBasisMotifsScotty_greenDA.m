% Define base directories and parameters
filePreprocessed = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed';
fileKeyword = '_green_dff_combined.mat'; % Suffix of processed imaging files
parameter_class = 'general_params_dual'; % Parameter file class used for local/remote bucket info
basis_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData/clusterW_output_DAmotifs_072125.mat';  % Path to basis motifs

% Search for all animal folders containing dual-imaging data
mListC = GrabFiles_sort_trials('m*jRGECO', 0, {'/Volumes/buschman/Rodent Data/dualImaging_parkj'});

% Iterate over each animal
for j = 1:numel(mListC)

    % Extract mouse ID (e.g., m1045) from folder name
    mIdC = regexp(mListC{j}, 'm\d{4}', 'match');
    if ~isempty(mIdC)
        mId = mIdC{1};
    else
        mId = '';
    end

    % Find session folders within the animal folder
    filePathSessions = find_keyword_containing_folder(mListC{j}, mId, 'recursive', false);

    % Iterate over each session
    for jj = 1:numel(filePathSessions)

        % Locate the _img folder for the session (used to define session identity)
        filePathImg = find_keyword_containing_folder(fullfile(filePathSessions{jj}, 'task'), '_img', 'recursive', false);
        if ~isempty(filePathImg)
            if iscell(filePathImg)
                filePathImg = filePathImg{1};
            end

            [~, header] = fileparts(filePathImg);  % Extract folder name as session identifier

            % Find corresponding preprocessed data folder using session name
            filePathPreprop = find_keyword_containing_folder(filePreprocessed, header, 'recursive', false);

            % Only proceed if preprocessed data exists
            if ~isempty(filePathPreprop)
                % Define where to save refit results
                save_dir = fullfile(fileparts(filePathImg), 'Matfiles');

                % Call function to submit chunked refitting jobs to Scotty
                Scotty_RefitBasisMotifs_Swarm_JP(filePathImg, fileKeyword, basis_dir, parameter_class, save_dir)
            end

            fprintf("Submitted jobs for session#%d of %s\n", jj, mId)
        end
    end
end
