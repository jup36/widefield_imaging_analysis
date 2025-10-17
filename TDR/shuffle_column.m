function xsh = shuffle_column(x, mode, trialIdx, blockSize)
% Shuffle a single column according to mode.
switch lower(string(mode))
    case "global"
        xsh = x(randperm(numel(x)));
    case "within_trial"
        assert(~isempty(trialIdx), 'RowIndex required for within_trial mode.');
        xsh = x;
        uTrials = unique(trialIdx(:))';
        for t = uTrials
            rows = find(trialIdx==t);
            xsh(rows) = x(rows(randperm(numel(rows))));
        end
    case "circular"
        % random circular shift
        sh = randi([0 numel(x)-1],1,1);
        xsh = circshift(x, sh);
    case "block"
        % shuffle blocks of length blockSize
        n = numel(x);
        idx = 1:blockSize:n;
        blocks = arrayfun(@(i) i:min(i+blockSize-1,n), idx,'uni',0);
        order = randperm(numel(blocks));
        xsh = zeros(n,1);
        pos = 1;
        for k = 1:numel(order)
            b = blocks{order(k)};
            m = numel(b);
            xsh(pos:pos+m-1) = x(b);
            pos = pos + m;
        end
    otherwise
        error('Unknown ShuffleMode: %s', mode);
end
end