function [featMat, sessInfo] = buildSessionFeaturesFromXcorr_dir(xcorrPosLagMatC, mListC, headerC)
% Build session-wise feature matrix from directional motif xcorr matrices.
%
% Each session's 30x30 matrix A is interpreted as A(i,j) = j after i
% (positive-lag correlation). We:
%   1) Fisher z-transform A
%   2) Drop the diagonal (i == j)
%   3) Flatten the remaining K*(K-1) entries into a feature vector
%
% OUTPUTS
%   featMat : [S x D] double, S = #sessions, D = K*(K-1)
%   sessInfo: table with session metadata

    nMice   = size(xcorrPosLagMatC,1);
    nSessMx = size(xcorrPosLagMatC,2);

    featMat_list    = {};
    mouseIdx_list   = [];
    mouseId_list    = {};
    sessWithin_list = [];
    header_list     = {};

    for i = 1:nMice
        thisPath = mListC{i};

        % Try to extract something like 'm1045' or 'm1893'
        tok = regexp(thisPath, '(m\d{3,5})', 'tokens', 'once');

        if isempty(tok)
            % Fallback: label by index if regex fails
            mId = sprintf('mouse_%02d', i);
            warning('buildSessionFeaturesFromXcorr_dir:NoMouseId', ...
                'Could not extract mouseId from "%s"; using "%s" instead.', ...
                thisPath, mId);
        else
            mId = tok{1};   % tok is a 1x1 cell containing the match
        end

        sessCount = 0;

        for j = 1:nSessMx
            A = xcorrPosLagMatC{i,j};
            if isempty(A)
                continue;
            end
            sessCount = sessCount + 1;

            % --- Fisher z-transform (clip to avoid Inf) ---
            A = max(min(A, 0.9999), -0.9999);
            Az = atanh(A);                % same size as A

            % --- drop diagonal, keep ALL off-diagonals (directional) ---
            K = size(Az,1);
            mask = ~eye(K);               % false on diag, true elsewhere

            featVec = Az(mask)';          % row vector, length K*(K-1)

            featMat_list{end+1,1} = featVec;      %#ok<AGROW>
            mouseIdx_list         = [mouseIdx_list;   i];         %#ok<AGROW>
            mouseId_list          = [mouseId_list;    {mId}];     %#ok<AGROW>
            sessWithin_list       = [sessWithin_list; sessCount]; %#ok<AGROW>
            header_list           = [header_list;     headerC{i,j}]; %#ok<AGROW>
        end
    end

    featMat = cell2mat(featMat_list);  % S x D

    sessInfo = table( ...
        (1:size(featMat,1))', ...
        mouseIdx_list, ...
        mouseId_list, ...
        sessWithin_list, ...
        header_list, ...
        'VariableNames', {'sessIdx','mouseIdx','mouseId','sessWithin','header'});
end
