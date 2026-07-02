function [tClean, yClean] = sanitizeTimeAndSignal(t, y)
% Make interp1-safe vectors:
%   - finite time and finite signal only
%   - sorted by time
%   - duplicate timestamps removed

t = double(t(:)');
y = double(y(:)');

n = min(numel(t), numel(y));
t = t(1:n);
y = y(1:n);

validI = isfinite(t) & isfinite(y);

tClean = t(validI);
yClean = y(validI);

if isempty(tClean)
    return;
end

[tClean, sortI] = sort(tClean, 'ascend');
yClean = yClean(sortI);

[tClean, uniqueI] = unique(tClean, 'stable');
yClean = yClean(uniqueI);

end
