% Helper to try add one continuous var (pads empties with NaNs first)
function [X_blocks, X_names] = try_add_cont(tbytDat, varField, outName, zThres, origCTS, X_blocks, X_names)
    if ~isfield(tbytDat, varField), return; end

    targetLen = numel(origCTS);  % length you expect per trial (e.g., -1:0.02:6)

    % Build a cell array of 1×targetLen rows; empty/missing → NaNs
    sigC = cell(1, numel(tbytDat));
    for ii = 1:numel(tbytDat)
        if ~isempty(tbytDat(ii).(varField))
            v = tbytDat(ii).(varField);
            v = v(1,:);                           % first row = signal
        else
            v = nan(1, targetLen);                % empty trial → NaNs
        end
        sigC{ii} = v;
    end

    % Convert to design column (keeps NaNs; timeseries_to_design can smooth/clip)
    Xc = timeseries_to_design(sigC, origCTS, ...
            'Epoch', [-0.9 5], 'Win', 0.1, 'Step', 0.05, ...
            'Method','mean', 'Fill', NaN, 'zscore', false, ... % <-- no z here
            'clipExtremes', true, 'zThres', zThres);

    Xc = Xc(:);  % column

    % Guard: append only if sizes match current design block
    if ~isempty(Xc) && size(Xc,1) == size(X_blocks{1},1)
        X_blocks{end+1} = Xc; %#ok<AGROW>
        X_names{end+1}  = outName; %#ok<AGROW>
    end
end