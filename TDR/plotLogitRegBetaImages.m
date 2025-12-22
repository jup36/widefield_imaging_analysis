function plotLogitRegBetaImages(logit_BetaC, glmRezC, glmLabelC, headerC)
% plotLogitRegBetaImages
%   Make and save images of time-resolved logistic regression beta weights
%   for one animal (multiple sessions).
%
% Inputs
%   logit_BetaC : 1 x nSess cell, each [nMotifs x T] beta matrix (decBins.W)
%   glmRezC     : 1 x nSess cell, each glmRez struct (needs .decBins.time)
%   glmLabelC   : 1 x nSess cell, each full path to glmRez file (for taskDir)
%   headerC     : 1 x nSess cell, each session header (e.g. m1045_122424)
%
% Junchol Park, 2025-07-31

for j = 1:numel(logit_BetaC)  % sessions
    if isempty(logit_BetaC{j}) || isempty(glmRezC{j}) || isempty(glmLabelC{j})
        continue
    end

    % ---- directories ----
    taskDirTok = regexp(glmLabelC{j}, '^(.*?task)[\\/]', 'tokens','once');
    if isempty(taskDirTok); continue; end
    taskDir = taskDirTok{1};

    figSaveDir_logitReg = fullfile(taskDir, "Figure", "logitRegGoNogo_H");
    if ~exist(figSaveDir_logitReg, "dir")
        mkdir(figSaveDir_logitReg);
    end

    % ---- data ----
    periCueTime = glmRezC{j}.decBins.time;   % 1 x T (or T x 1)
    periCueTime = periCueTime(:)';           % force row
    B           = logit_BetaC{j};            % nMotifs x T
    nMotifs     = size(B,1);

    % ---- plot ----
    hFig_logic = figure('Visible','off');  %#ok<NASGU>  % off = faster when batching

    imagesc(periCueTime, 1:nMotifs, B);
    title(['logitRegGoNogoH_', headerC{j}], 'Interpreter','none');
    clim([-.2 .2]);
    ax = gca;
    set(ax, 'YDir','normal', ...    % motif 1 at bottom
            'TickDir','out');

    xlabel('Time from cue (s)');
    ylabel('Motif');

    % ---- x-ticks: -0.5 : 0.5 : 5 (only those within data range) ----
    desiredTicks = -0.5:0.5:5;
    validTicks   = desiredTicks(desiredTicks >= min(periCueTime) & ...
                                desiredTicks <= max(periCueTime));
    xticks(validTicks);
    xticklabels(compose('%.1f', validTicks));

    % ---- vertical white dotted lines at 0, 2, 4 s (if in range) ----
    hold on
    lineTimes = [0 2 4];
    for lt = lineTimes
        if lt >= min(periCueTime) && lt <= max(periCueTime)
            [~, idx] = min(abs(periCueTime - lt));
            xline(periCueTime(idx), 'w:', 'LineWidth', 1.5);
        end
    end
    hold off

    % ---- save ----
    outName = fullfile(figSaveDir_logitReg, ...
        ['logitRegGoNogoH_', headerC{j}]);
    print(gcf, outName, '-dpdf', '-painters', '-bestfit');

    close(gcf);
end
end
