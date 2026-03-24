function stats = runGroupTimeRMANOVA(group1, group2, varargin)
%RUNGROUPTIMERMANOVA Repeated-measures ANOVA for group x time data.
%
% stats = runGroupTimeRMANOVA(group1, group2, 'Name', value, ...)
%
% INPUT
%   group1 : [n1 x nTime] numeric matrix
%            rows = subjects/animals, columns = time bins
%   group2 : [n2 x nTime] numeric matrix
%            rows = subjects/animals, columns = time bins
%
% NAME-VALUE PAIRS
%   'groupNames'        : 1x2 cell/string array of group names
%                         default = {'group1','group2'}
%   'timeX'             : 1 x nTime vector of time values
%                         default = 1:nTime
%   'timeBounds'        : [] (default) or [tMin tMax]
%                         If provided, only bins with timeX in [tMin, tMax]
%                         are retained for all analyses.
%   'alpha'             : significance threshold
%                         default = 0.05
%   'multcompareMethod' : 'fdr' | 'holm'
%                         default = 'fdr'
%   'doPlot'            : logical scalar
%                         default = true
%   'plotTitle'         : char/string
%                         default = ''
%
% OUTPUT
%   stats : struct containing results

%% Parse inputs
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'group1', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(ip, 'group2', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));

addParameter(ip, 'groupNames', {'group1','group2'}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'timeX', [], @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'timeBounds', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
addParameter(ip, 'alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(ip, 'multcompareMethod', 'fdr', @(x) ischar(x) || isstring(x));
addParameter(ip, 'doPlot', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'plotTitle', '', @(x) ischar(x) || isstring(x));

parse(ip, group1, group2, varargin{:});
P = ip.Results;

groupNames = cellstr(string(P.groupNames));
if numel(groupNames) ~= 2
    error('groupNames must contain exactly two entries.');
end

[n1, nTime] = size(group1);
[n2, nTime2] = size(group2);
if nTime ~= nTime2
    error('group1 and group2 must have the same number of columns (time bins).');
end

if isempty(P.timeX)
    timeX = 1:nTime;
else
    timeX = P.timeX(:)';
    if numel(timeX) ~= nTime
        error('numel(timeX) must match number of time bins.');
    end
end

method = lower(string(P.multcompareMethod));
if ~ismember(method, ["fdr","holm"])
    error('multcompareMethod must be ''fdr'' or ''holm''.');
end

%% Apply temporal indexing before all statistics
timeMask = true(1, numel(timeX));

if ~isempty(P.timeBounds)
    tMin = P.timeBounds(1);
    tMax = P.timeBounds(2);
    timeMask = (timeX >= tMin) & (timeX <= tMax);

    if ~any(timeMask)
        error('timeBounds excluded all time bins.');
    end

    group1 = group1(:, timeMask);
    group2 = group2(:, timeMask);
    timeX  = timeX(timeMask);
end

% Update nTime after temporal selection
nTime = size(group1, 2);

%% Build repeated-measures table
Y = double([group1; group2]);
nSubj = size(Y,1);

group = [repmat(string(groupNames{1}), n1, 1); ...
         repmat(string(groupNames{2}), n2, 1)];
subject = (1:nSubj)';

dataTbl = array2table(Y, 'VariableNames', cellstr("t" + string(1:nTime)));
T = table(subject, categorical(group), 'VariableNames', {'Subj','Group'});
T = [T dataTbl];

within = table((1:nTime)', 'VariableNames', {'Time'});

%% Fit repeated-measures model
rm = fitrm(T, sprintf('t1-t%d ~ Group', nTime), 'WithinDesign', within);

%% Repeated-measures ANOVA
ranovaTbl = ranova(rm, 'WithinModel', 'Time');

% Keep for inspection only
try
    betweenTbl = anova(rm);
catch
    betweenTbl = table();
end

%% Extract Time and Group x Time from ranova
[~, pTime, pGroupTime, effectInfo] = local_extract_effect_ps(ranovaTbl, betweenTbl);

%% Compute Group main effect directly from subject-wise means across time
subjMean1 = mean(group1, 2, 'omitnan');
subjMean2 = mean(group2, 2, 'omitnan');

[~, pGroup, ~, groupStats] = ttest2(subjMean1, subjMean2, 'Vartype', 'unequal');

%% Per-bin post hoc group comparisons
p_unc = nan(1, nTime);
tstat = nan(1, nTime);
df    = nan(1, nTime);

for t = 1:nTime
    x1 = double(group1(:,t));
    x2 = double(group2(:,t));

    [~, p, ~, s] = ttest2(x1, x2, 'Vartype', 'unequal');

    p_unc(t) = p;
    tstat(t) = s.tstat;
    df(t)    = s.df;
end

%% Multiple-comparison correction
switch method
    case "fdr"
        [p_corr, sig] = local_fdr_bh(p_unc, P.alpha);
    case "holm"
        [p_corr, sig] = local_holm_bonferroni(p_unc, P.alpha);
end

%% Collect output
stats = struct();
stats.rm = rm;
stats.ranovaTbl = ranovaTbl;
stats.betweenTbl = betweenTbl;

stats.pGroup = pGroup;
stats.pTime = pTime;
stats.pGroupTime = pGroupTime;
stats.effectInfo = effectInfo;

stats.subjMean1 = subjMean1;
stats.subjMean2 = subjMean2;
stats.groupStats = groupStats;

stats.p_unc = p_unc;
stats.p_corr = p_corr;
stats.sig = sig;
stats.tstat = tstat;
stats.df = df;

stats.timeX = timeX;
stats.timeBounds = P.timeBounds;
stats.timeMask = timeMask;

stats.groupNames = groupNames;
stats.alpha = P.alpha;
stats.multcompareMethod = char(method);

stats.nGroup1 = n1;
stats.nGroup2 = n2;
stats.meanGroup1 = mean(group1, 1, 'omitnan');
stats.meanGroup2 = mean(group2, 1, 'omitnan');
stats.semGroup1 = std(group1, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(group1),1));
stats.semGroup2 = std(group2, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(group2),1));

%% Optional display
fprintf('\n=== Repeated-measures ANOVA ===\n');
if ~isempty(P.timeBounds)
    fprintf('Time bounds:         [%.4g %.4g]\n', P.timeBounds(1), P.timeBounds(2));
    fprintf('Retained bins:       %d\n', nTime);
end
fprintf('Group effect:        p = %.4g\n', stats.pGroup);
fprintf('Time effect:         p = %.4g\n', stats.pTime);
fprintf('Group x Time effect: p = %.4g\n', stats.pGroupTime);
fprintf('Per-bin correction:  %s\n', stats.multcompareMethod);
fprintf('Significant bins:    %d / %d\n', sum(stats.sig), nTime);

%% Optional plotting
if P.doPlot
    local_plot_results(stats, P.plotTitle);
end

end


%% ---------- Local helper functions ----------

function [pGroup, pTime, pGroupTime, info] = local_extract_effect_ps(ranovaTbl, betweenTbl)

pGroup = NaN;
pTime = NaN;
pGroupTime = NaN;

info = struct();
info.ranova_pcol = '';
info.between_pcol = '';
info.ranova_rows = [];
info.between_rows = [];
info.ranova_varnames = string(ranovaTbl.Properties.VariableNames);

if ~isempty(betweenTbl)
    info.between_varnames = string(betweenTbl.Properties.VariableNames);
else
    info.between_varnames = string.empty;
end

%% ----- ranova table: Time and Group x Time -----
pColR = '';
if any(strcmpi(ranovaTbl.Properties.VariableNames, 'pValue'))
    pColR = ranovaTbl.Properties.VariableNames{strcmpi(ranovaTbl.Properties.VariableNames, 'pValue')};
elseif any(strcmpi(ranovaTbl.Properties.VariableNames, 'pValueGG'))
    pColR = ranovaTbl.Properties.VariableNames{strcmpi(ranovaTbl.Properties.VariableNames, 'pValueGG')};
elseif any(strcmpi(ranovaTbl.Properties.VariableNames, 'pValueHF'))
    pColR = ranovaTbl.Properties.VariableNames{strcmpi(ranovaTbl.Properties.VariableNames, 'pValueHF')};
elseif any(strcmpi(ranovaTbl.Properties.VariableNames, 'pValueLB'))
    pColR = ranovaTbl.Properties.VariableNames{strcmpi(ranovaTbl.Properties.VariableNames, 'pValueLB')};
end
info.ranova_pcol = pColR;

rowNamesR = string(ranovaTbl.Properties.RowNames);
rowNamesR_clean = local_clean_terms(rowNamesR);
info.ranova_rows = table(rowNamesR(:), rowNamesR_clean(:), ...
    'VariableNames', {'Original','Cleaned'});

if ~isempty(pColR)
    iTime = find(contains(rowNamesR_clean, 'time') & ...
                 ~contains(rowNamesR_clean, 'group:time') & ...
                 ~contains(rowNamesR_clean, 'time:group') & ...
                 ~contains(rowNamesR_clean, 'errortime') & ...
                 ~contains(rowNamesR_clean, 'time:error'), 1, 'first');

    iGroupTime = find(contains(rowNamesR_clean, 'group:time') | ...
                      contains(rowNamesR_clean, 'time:group') | ...
                      contains(rowNamesR_clean, 'group*time') | ...
                      contains(rowNamesR_clean, 'time*group'), 1, 'first');

    if ~isempty(iTime)
        pTime = ranovaTbl.(pColR)(iTime);
    end
    if ~isempty(iGroupTime)
        pGroupTime = ranovaTbl.(pColR)(iGroupTime);
    end
end

end


function cleaned = local_clean_terms(x)
x = string(x);
cleaned = lower(strtrim(x));
cleaned = replace(cleaned, " ", "");
cleaned = replace(cleaned, "*", ":");
cleaned = replace(cleaned, "×", ":");
cleaned = replace(cleaned, "(between)", "");
cleaned = replace(cleaned, "(within)", "");
cleaned = replace(cleaned, "(", "");
cleaned = replace(cleaned, ")", "");
end


function [p_adj, sig] = local_fdr_bh(p, alpha)
p = p(:);
m = numel(p);

[ps, idx] = sort(p);
thresh = (1:m)'/m * alpha;

below = ps <= thresh;
sig = false(m,1);

if any(below)
    k = find(below, 1, 'last');
    sig(idx(1:k)) = true;
end

p_adj_sorted = nan(m,1);
for i = 1:m
    p_adj_sorted(i) = min(1, ps(i) * m / i);
end
for i = m-1:-1:1
    p_adj_sorted(i) = min(p_adj_sorted(i), p_adj_sorted(i+1));
end

p_adj = nan(m,1);
p_adj(idx) = p_adj_sorted;

p_adj = p_adj(:)';
sig = sig(:)';
end


function [p_adj, sig] = local_holm_bonferroni(p, alpha)
p = p(:);
m = numel(p);

[ps, idx] = sort(p);
sig_sorted = false(m,1);

for i = 1:m
    if ps(i) <= alpha / (m - i + 1)
        sig_sorted(i) = true;
    else
        break;
    end
end

sig = false(m,1);
sig(idx(sig_sorted)) = true;

p_adj_sorted = nan(m,1);
for i = 1:m
    p_adj_sorted(i) = min(1, (m - i + 1) * ps(i));
end
for i = 2:m
    p_adj_sorted(i) = max(p_adj_sorted(i), p_adj_sorted(i-1));
end

p_adj = nan(m,1);
p_adj(idx) = p_adj_sorted;

p_adj = p_adj(:)';
sig = sig(:)';
end


function local_plot_results(stats, plotTitle)
figure('Color', 'w');

subplot(3,1,1);
plot(stats.timeX, stats.meanGroup1, 'LineWidth', 1.5); hold on;
plot(stats.timeX, stats.meanGroup2, 'LineWidth', 1.5);
xline(0, '--k');
legend(stats.groupNames, 'Location', 'best', 'Interpreter', 'none');
ylabel('Mean');

titleStr = string(plotTitle);
if strlength(titleStr) == 0
    title(sprintf('Group p=%.3g, Time p=%.3g, Group×Time p=%.3g', ...
        stats.pGroup, stats.pTime, stats.pGroupTime));
else
    title(sprintf('%s | Group p=%.3g, Time p=%.3g, Group×Time p=%.3g', ...
        titleStr, stats.pGroup, stats.pTime, stats.pGroupTime), ...
        'Interpreter', 'none');
end

subplot(3,1,2);
plot(stats.timeX, stats.p_corr, 'k', 'LineWidth', 1);
hold on;
yline(stats.alpha, '--r');
xline(0, '--k');
ylabel('Corrected p');
title(sprintf('Per-bin post hoc (%s corrected)', stats.multcompareMethod), ...
    'Interpreter', 'none');

subplot(3,1,3);
stem(stats.timeX(stats.sig), ones(1,sum(stats.sig)), 'filled');
hold on;
xline(0, '--k');
ylim([0 1.2]);
xlabel('Time');
ylabel('sig');
title('Significant bins');
end