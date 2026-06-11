function writtenFiles = write_lk_parameter_class(baseParamClass, newParamClass, L, K)
% Write a new parameter class file based on an existing parameter class.
%
% Writes the generated class file to:
%   1) every folder where baseParamClass exists on the MATLAB path
%   2) the Slurm-accessible ParameterClasses folder
%
% This handles the case where the local copy is shadowed by the Z: copy.

    %% Locate all visible copies of base parameter class
    baseFiles = which(baseParamClass, '-all');

    if isempty(baseFiles)
        error('Could not find base parameter class file on MATLAB path: %s', baseParamClass);
    end

    if ischar(baseFiles) || isstring(baseFiles)
        baseFiles = cellstr(baseFiles);
    end

    % Clean possible display annotations after " % "
    % This is defensive; usually which(...,'-all') returns clean paths.
    for ii = 1:numel(baseFiles)
        thisFile = baseFiles{ii};
        pctIdx = strfind(thisFile, ' % ');
        if ~isempty(pctIdx)
            thisFile = strtrim(thisFile(1:pctIdx(1)-1));
        end
        baseFiles{ii} = thisFile;
    end

    fprintf('\nAll detected copies of %s:\n', baseParamClass);
    for ii = 1:numel(baseFiles)
        fprintf('  %s\n', baseFiles{ii});
    end

    % Use the first copy as the template source.
    % This matches MATLAB's active class resolution.
    templateFile = baseFiles{1};

    %% Get all folders containing baseParamClass
    baseParamDirs = cell(size(baseFiles));

    for ii = 1:numel(baseFiles)
        [baseParamDirs{ii}, ~, ~] = fileparts(baseFiles{ii});
    end

    %% Fixed Slurm-accessible ParameterClasses directory
    slurmParamDir = compatiblepath( ...
        'Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\ParameterClasses');

    %% Destination folders: all base-class folders + Slurm folder
    destDirs = [baseParamDirs(:); {slurmParamDir}];

    % Remove duplicate folders while preserving order
    destDirsNorm = cellfun(@(x) lower(char(x)), destDirs, 'UniformOutput', false);
    [~, uniqueIdx] = unique(destDirsNorm, 'stable');
    destDirs = destDirs(uniqueIdx);

    fprintf('\nGenerated parameter class will be written to:\n');
    for dd = 1:numel(destDirs)
        fprintf('  %s\n', destDirs{dd});
    end

    %% Check all destination folders before writing anything
    for dd = 1:numel(destDirs)

        thisDir = destDirs{dd};

        if exist(thisDir, 'dir') ~= 7
            error(['Required parameter-class output folder is not accessible:\n%s\n\n' ...
                   'Cannot write generated parameter class.'], thisDir);
        end

        testFile = fullfile(thisDir, sprintf('__write_test_%s.tmp', ...
            datestr(now, 'yyyymmdd_HHMMSSFFF')));

        fid = fopen(testFile, 'w');

        if fid == -1
            error(['Parameter-class output folder exists but is not writable:\n%s\n\n' ...
                   'Check drive mounting, permissions, or network access.'], thisDir);
        end

        fprintf(fid, 'write test\n');
        fclose(fid);

        if exist(testFile, 'file') ~= 2
            error('Write test failed unexpectedly in folder:\n%s', thisDir);
        end

        delete(testFile);

    end

    %% Read template class file
    txt = fileread(templateFile);

    %% Replace classdef name
    txt = regexprep( ...
        txt, ...
        ['classdef\s+', baseParamClass], ...
        ['classdef ', newParamClass], ...
        'once');

    %% Replace K assignment
    txt = regexprep( ...
        txt, ...
        '(\n\s*K\s*=\s*)[0-9]+(\s*[;\n])', ...
        sprintf('$1%d$2', K), ...
        'once');

    %% Replace L assignment
    txt = regexprep( ...
        txt, ...
        '(\n\s*L\s*=\s*)[0-9]+(\s*[;\n])', ...
        sprintf('$1%d$2', L), ...
        'once');

    %% Safety checks
    if ~contains(txt, ['classdef ', newParamClass])
        error('Failed to update classdef name from %s to %s.', ...
            baseParamClass, newParamClass);
    end

    if ~contains(txt, sprintf('K = %d', K)) && ~contains(txt, sprintf('K=%d', K))
        warning('Could not verify K replacement for generated class: %s', newParamClass);
    end

    if ~contains(txt, sprintf('L = %d', L)) && ~contains(txt, sprintf('L=%d', L))
        warning('Could not verify L replacement for generated class: %s', newParamClass);
    end

    %% Write generated class to all destination folders
    writtenFiles = cell(numel(destDirs), 1);

    for dd = 1:numel(destDirs)

        thisFile = fullfile(destDirs{dd}, [newParamClass, '.m']);

        fid = fopen(thisFile, 'w');

        if fid == -1
            error('Could not open generated parameter class for writing:\n%s', thisFile);
        end

        fwrite(fid, txt);
        fclose(fid);

        if exist(thisFile, 'file') ~= 2
            error('File write appeared to succeed, but file was not found afterward:\n%s', thisFile);
        end

        writtenFiles{dd} = thisFile;

    end

    %% Refresh MATLAB path/class cache
    for dd = 1:numel(destDirs)
        addpath(destDirs{dd});
    end

    rehash;

    %% Report
    fprintf('\nParameter class written successfully to:\n');
    for dd = 1:numel(writtenFiles)
        fprintf('  %s\n', writtenFiles{dd});
    end

end