function inspect_cp_components(cpRez, sessionInfo, timeVec, R)
% INSPECT_CP_COMPONENTS  Quick plots for a chosen CP rank.
%
%   inspect_cp_components(cpRez, sessionInfo, timeVec, R)
%
% INPUT
%   cpRez       : struct returned by run_cp_rank
%   sessionInfo : table with fields
%                   .sessIdx, .mouseIdx, .sessWithinMouse, .mouseId, ...
%   timeVec     : [T x 1] or [1 x T] time axis for decoder weights
%   R           : rank to visualize (must be in cpRez.R)
%
% OUTPUT
%   (figures only; no returns)
%
% Plots:
%   1) Session embedding (comp 1 vs comp 2), colored by mouse.
%   2) Session trajectories along comp 1 vs session # within mouse.
%   3) Motif loadings per component.
%   4) Temporal profile per component.

    % ---- find index for requested rank ----
    idxR = find(cpRez.R == R, 1);
    if isempty(idxR)
        error('Rank R = %d not found in cpRez.R.', R);
    end

    U = cpRez.factors{idxR};
    U_sess  = U{1};   % S x R
    M_motif = U{2};   % K x R
    T_time  = U{3};   % T x R

    [S, Rchk] = size(U_sess);
    if Rchk ~= R
        warning('U_sess size mismatch (columns=%d, R=%d).', Rchk, R);
    end

    if numel(timeVec) ~= size(T_time,1)
        warning('timeVec length (%d) != #rows in T_time (%d).', ...
            numel(timeVec), size(T_time,1));
    end
    timeVec = timeVec(:)';  % row

    mouseIdsUnique = unique(sessionInfo.mouseIdx);
    cmap = lines(numel(mouseIdsUnique));

    %% 1) Session embedding (component 1 vs 2)
    h1 = figure('Color','w'); hold on;
    for ii = 1:numel(mouseIdsUnique)
        m   = mouseIdsUnique(ii);
        msk = (sessionInfo.mouseIdx == m);
        scatter(U_sess(msk,1), U_sess(msk,2), 60, cmap(ii,:), ...
            'filled', 'DisplayName', sessionInfo.mouseId{find(msk,1)});
        plot(U_sess(msk,1), U_sess(msk,2), '-', 'Color', cmap(ii,:));
    end
    xlabel('CP component 1'); ylabel('CP component 2');
    title(sprintf('Session embedding (R = %d)', R));
    legend('Location','bestoutside'); box on; grid on;

    %% 2) Session trajectories along component 1
    h2 = figure('Color','w'); hold on;
    for ii = 1:numel(mouseIdsUnique)
        m   = mouseIdsUnique(ii);
        msk = (sessionInfo.mouseIdx == m);
        % sort by within-mouse session index
        [~,ord] = sort(sessionInfo.sessWithinMouse(msk));
        idxSess_m = find(msk); idxSess_m = idxSess_m(ord);

        plot(sessionInfo.sessWithinMouse(idxSess_m), ...
             U_sess(idxSess_m,1), '-o', ...
             'Color', cmap(ii,:), ...
             'DisplayName', sessionInfo.mouseId{idxSess_m(1)});
    end
    xlabel('Session # within mouse');
    ylabel('CP component 1 score');
    title(sprintf('Learning trajectories along component 1 (R = %d)', R));
    legend('Location','bestoutside'); box on; grid on;

    %% 3) Motif loadings
    [K, ~] = size(M_motif);
    h3 = figure('Color','w');
    for r = 1:R
        subplot(R,1,r);
        stem(1:K, M_motif(:,r), 'filled');
        xlabel('Motif #'); ylabel(sprintf('Comp %d', r));
        title(sprintf('Motif factor for component %d', r));
        box off; grid on;
    end

    %% 4) Temporal profiles
    h4 = figure('Color','w'); hold on;
    for r = 1:R
        plot(timeVec, T_time(:,r), 'LineWidth', 1.5);
    end
    xlabel('Time (s)');
    ylabel('Temporal factor');
    title(sprintf('Temporal profiles of components (R = %d)', R));
    legend(arrayfun(@(r)sprintf('Comp %d',r),1:R,'UniformOutput',false), ...
           'Location','best'); 
    grid on; box on;

    % return handles if you want to capture them:
    if nargout > 0
        varargout{1} = struct('hSessEmbed',h1,'hTraj',h2,'hMotif',h3,'hTime',h4);
    end
end
