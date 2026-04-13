function plotDataCHistogram(dataC, nbins)

% plotDataCHistogram(dataC, nbins)
% dataC : cell array (e.g., 15x2), each cell = 3D stack (64x64xT)
% nbins : number of histogram bins (default = 200)

if nargin < 2
    nbins = 200;
end

allVals = [];  % accumulator

for i = 1:numel(dataC)
    if isempty(dataC{i})
        continue;
    end
    
    dat = dataC{i};
    
    % Flatten and remove NaNs
    dat = dat(:);
    dat = dat(~isnan(dat));
    
    allVals = [allVals; dat]; %#ok<AGROW>
end

% ---- Plot histogram ----
figure;
histogram(allVals, nbins);
xlabel('Value');
ylabel('Count');
title('Distribution of all non-NaN values in dataC');

% ---- Print summary stats ----
fprintf('--- Summary stats ---\n');
fprintf('Min: %.4f\n', min(allVals));
fprintf('Max: %.4f\n', max(allVals));
fprintf('Mean: %.4f\n', mean(allVals));
fprintf('Std: %.4f\n', std(allVals));

% ---- Show percentiles (very useful for outliers) ----
prc = prctile(allVals, [0.1 1 5 50 95 99 99.9]);
fprintf('Percentiles [0.1 1 5 50 95 99 99.9]:\n');
disp(prc);

end