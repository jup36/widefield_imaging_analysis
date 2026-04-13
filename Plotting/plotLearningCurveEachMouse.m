function hC = plotLearningCurveEachMouse(dat, cMat, figSaveDir, fileTag)
% plotLearningCurveEachMouse
% Plot one learning-curve figure per mouse.
%
% Inputs:
%   dat        - Nx2 cell array
%                dat{i,1}: 1xD cell array of headers (e.g., {'m1045_122424', ...})
%                dat{i,2}: 1xD numeric array of values
%   cMat       - Nx3 colormap matrix, one RGB row per mouse
%   figSaveDir - optional directory for saving figures
%   fileTag    - optional string appended to output filename
%
% Output:
%   hC         - Nx1 cell array of figure handles
%
% Example:
%   hC = plotLearningCurveEachMouse(dPrmTcollectC, cMat, figSaveDir, 'dPrm');

    if nargin < 2 || isempty(cMat)
        cMat = lines(size(dat,1));
    end
    if nargin < 3
        figSaveDir = [];
    end
    if nargin < 4 || isempty(fileTag)
        fileTag = 'learningCurve';
    end

    numMice = size(dat, 1);
    hC = cell(numMice, 1);

    for iMouse = 1:numMice
        mouseHeader = dat{iMouse, 1};
        mouseData   = dat{iMouse, 2};

        if isempty(mouseHeader) || isempty(mouseData)
            warning('Mouse index %d has empty header or data. Skipping.', iMouse);
            continue;
        end

        % Ensure row format
        if iscolumn(mouseData)
            mouseData = mouseData';
        end
        if iscolumn(mouseHeader)
            mouseHeader = mouseHeader';
        end

        numDays = length(mouseData);
        xDays   = 1:numDays;

        % Extract mouse ID from first header
        mId = regexp(mouseHeader{1}, 'm\d{4,5}', 'match', 'once');
        if isempty(mId)
            mId = sprintf('mouse%02d', iMouse);
        end

        % Make figure
        h = figure('Color', 'w');
        hold on;

        plot(xDays, mouseData, '-o', ...
            'Color', cMat(iMouse, :), ...
            'MarkerFaceColor', cMat(iMouse, :), ...
            'MarkerEdgeColor', 'none', ...
            'LineWidth', 1.5);

        % x-axis labels = headers, with literal interpretation
        set(gca, 'XTick', xDays, ...
                 'XTickLabel', mouseHeader, ...
                 'XTickLabelRotation', 45, ...
                 'TickLabelInterpreter', 'none');

        xlim([0.5, numDays + 0.5]);

        xlabel('Session');
        ylabel('Performance');
        title(sprintf('Learning Curve: %s', mId), 'Interpreter', 'none');

        box off;
        grid on;
        set(gca, 'TickDir', 'out');

        hold off;

        hC{iMouse} = h;

        % Save if requested
        if ~isempty(figSaveDir)
            if ~exist(figSaveDir, 'dir')
                mkdir(figSaveDir);
            end

            saveName = fullfile(figSaveDir, sprintf('%s_%s', fileTag, mId));
            print(h, saveName, '-dpdf', '-vector', '-bestfit');
        end
    end
end