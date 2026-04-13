function dataC_clean = clipDataC_percentile(dataC, lowP, highP)

if nargin < 2, lowP = 0.1; end
if nargin < 3, highP = 99.9; end

% collect all values
allVals = [];
for i = 1:numel(dataC)
    if isempty(dataC{i}), continue; end
    v = dataC{i}(:);
    v = v(~isnan(v));
    allVals = [allVals; v]; %#ok<AGROW>
end

lo = prctile(allVals, lowP);
hi = prctile(allVals, highP);

fprintf('Clipping range: [%.3f, %.3f]\n', lo, hi);

% apply clipping
dataC_clean = dataC;
for i = 1:numel(dataC)
    if isempty(dataC{i}), continue; end
    dat = dataC{i};
    dat(dat < lo) = lo;
    dat(dat > hi) = hi;
    dataC_clean{i} = dat;
end

end