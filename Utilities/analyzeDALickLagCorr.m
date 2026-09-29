function out = analyzeDALickLagCorr(rezC, sessT, varargin)
%ANALYZEDALICKLAGCORR
%   Lag-resolved, trial-by-trial correlation between binned lick count and
%   binned DA, within session. Companion to analyzeDAVigorControl, which
%   collapsed each trial to one scalar per window; this keeps the time
%   structure and asks WHEN licking and DA covary.
%
%   THE STATISTIC
%   -------------
%   For a given DA bin b and lag L, take across trials of one session:
%       x = lickCount(b - L, :)      (lick L bins EARLIER than the DA bin)
%       y = daBinned(sig, b, :)
%   and correlate them ACROSS TRIALS. Positive lag therefore means LICKING
%   LEADS DA. Averaging Fisher-z over bins b within a window gives one r per
%   lag per session; sessions are then averaged within animal, and animals
%   tested across with the exact sign-flip test -- the same two-stage
%   structure used elsewhere in this project.
%
%   WHY ACROSS TRIALS, NOT ALONG TIME
%   ---------------------------------
%   Correlating the two time courses WITHIN a trial mostly recovers their
%   shared task locking: both rise after the tone, so any two such signals
%   correlate strongly whether or not they are related trial to trial. The
%   across-trial correlation at fixed bin asks the question that matters --
%   on trials where the animal licked more than usual at this moment, was DA
%   also higher than usual? -- and is what a vigor account actually predicts.
%
%   WHAT THE LAG PROFILE DISTINGUISHES
%   ----------------------------------
%   A movement/vigor contribution should peak at or just before zero lag,
%   with a width set by the GRAB-DA indicator (hundreds of ms), because DA
%   would be reporting the licking itself. A signal related to outcome or
%   learning need not peak there, and can lead the licking. The profile is
%   therefore more informative than any single correlation coefficient, and
%   it is the part the scalar analysis threw away.
%
%   NOTE ON WHAT THIS CANNOT SHOW: a positive peak near zero lag is
%   CONSISTENT with DA reporting licking, but licking and DA may both be
%   driven by the same trial-to-trial variable (arousal, expected value), so
%   it is not evidence that licking causes the DA. Combined with the
%   analyzeDAVigorControl result -- session-order growth surviving the lick
%   control -- a nonzero within-session correlation would mean the two
%   effects live in different variance, not that one explains the other.
%
%   out = analyzeDALickLagCorr(rezC, sessT, ...)
%
% NAME-VALUE
%   'signal'       : 'global' (default) or motif number 1..K
%   'trialTypes'   : {'hit','cr'}
%   'window'       : [2 4] s -- DA bins whose correlations are averaged
%   'maxLagBins'   : 6  -- lags spanned, in coarse bins (+/-)
%   'minTrials'    : 20 -- a session needs this many trials of that type
%   'byLearnBin'   : false -- true splits sessions into early/intermediate/
%                    late (needs sessT.learnBin) and returns one profile each
%   'fastLearners' / 'slowLearners' : cohort lists
%   'figSaveDir'   : ''
%
% OUTPUT
%   .perSession  r [nLag x 1] per session x trialType (Fisher-z averaged over bins)
%   .perAnimal   r [nLag x 1] per animal x trialType (mean over its sessions)
%   .stats       per trialType x lag: mean r across animals, sign-flip p
%   .lagSec      lag axis in seconds
%   .fig

%% ---- parse ----
p = inputParser;
p.addParameter('signal', 'global', @(x) (ischar(x)||isstring(x)) || (isnumeric(x)&&isscalar(x)));
p.addParameter('trialTypes', {'hit','cr'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('window', [2 4], @(x) isnumeric(x) && numel(x)==2 && x(2)>x(1));
p.addParameter('maxLagBins', 6, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('minTrials', 20, @(x) isnumeric(x) && isscalar(x) && x>=5);
p.addParameter('byLearnBin', false, @(x) islogical(x) && isscalar(x));
p.addParameter('fastLearners', {'m1044','m1045','m1092','m1094'}, @iscell);
p.addParameter('slowLearners', {'m1048','m1049','m1613','m1859','m1873'}, @iscell);
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

trialTypes = cellstr(opt.trialTypes);
maxL = opt.maxLagBins;
lags = -maxL:maxL;

if isnumeric(opt.signal)
    sigRow = opt.signal + 1;                        % row 1 is global
    sigLabel = sprintf('motif %d DA', opt.signal); sigTag = sprintf('motif%02d', opt.signal);
else
    assert(strcmpi(opt.signal,'global'), 'signal must be ''global'' or a motif number.');
    sigRow = 1; sigLabel = 'global DA'; sigTag = 'global';
end

useBin = opt.byLearnBin && ismember('learnBin', sessT.Properties.VariableNames);
if opt.byLearnBin && ~useBin
    warning('LagCorr:NoLearnBin', 'sessT has no learnBin column -- ignoring byLearnBin.');
end
binNames = {'early','intermediate','late'};

%% ---- per session ----
rows = struct('animal',{},'group',{},'header',{},'learnBin',{},'trialType',{}, ...
    'r',{},'nTrials',{},'nBins',{});
binCtr = []; dtBin = NaN;

for i = 1:height(sessT)
    aId = char(sessT.animal{i});
    if ismember(aId, opt.fastLearners), g = 'fast';
    elseif ismember(aId, opt.slowLearners), g = 'slow';
    else, continue;
    end
    r = rezC{sessT.a(i), sessT.s(i)};
    if ~isfield(r, 'trial'), continue; end

    if isempty(binCtr)
        binCtr = r.trial.binCenters;
        dtBin  = r.trial.coarseBinSec;
    end
    bSel = find(binCtr >= opt.window(1) & binCtr < opt.window(2));
    if isempty(bSel), continue; end

    lb = "";
    if useBin, lb = string(sessT.learnBin(i)); end

    for t = 1:numel(trialTypes)
        tt = trialTypes{t};
        m = r.trial.type == tt;
        n = sum(m);
        if n < opt.minTrials, continue; end

        LC = double(r.trial.lickCount(:, m));           % [B x n]
        DA = double(squeeze(r.trial.daBinned(sigRow, :, m)));   % [B x n]
        B = size(LC, 1);

        z = nan(numel(lags), 1);
        for li = 1:numel(lags)
            L = lags(li);
            acc = []; 
            for b = bSel(:)'
                bl = b - L;                              % lick bin, L bins earlier
                if bl < 1 || bl > B, continue; end
                x = LC(bl, :)'; y = DA(b, :)';
                ok = isfinite(x) & isfinite(y);
                if sum(ok) < opt.minTrials || std(x(ok)) == 0 || std(y(ok)) == 0, continue; end
                rr = corr(x(ok), y(ok));
                acc(end+1) = atanh(max(min(rr, 0.999), -0.999)); %#ok<AGROW>
            end
            if ~isempty(acc), z(li) = mean(acc); end
        end

        rows(end+1) = struct('animal', aId, 'group', g, 'header', char(sessT.header{i}), ...
            'learnBin', lb, 'trialType', tt, 'r', tanh(z)', 'nTrials', n, ...
            'nBins', numel(bSel)); %#ok<AGROW>
    end
end
assert(~isempty(rows), 'No usable sessions -- has the batch been rerun so rezC has r.trial?');
perSession = struct2table(rows);
lagSec = lags * dtBin;

%% ---- per animal, then across animals ----
groupsToRun = {''};
if useBin, groupsToRun = binNames; end

perAnimal = struct('animal',{},'group',{},'trialType',{},'learnBin',{},'r',{},'nSessions',{});
stats = table();

for t = 1:numel(trialTypes)
    tt = trialTypes{t};
    for gI = 1:numel(groupsToRun)
        lb = groupsToRun{gI};
        sel = strcmp(perSession.trialType, tt);
        if useBin, sel = sel & perSession.learnBin == lb; end
        Ssel = perSession(sel, :);
        if isempty(Ssel), continue; end

        aList = unique(Ssel.animal, 'stable');
        A = nan(numel(aList), numel(lags));
        for a = 1:numel(aList)
            ra = rowsOf(Ssel.r, strcmp(Ssel.animal, aList{a}));   % [nSess x nLag]
            % Fisher-z average over that animal's sessions, equally weighted
            A(a, :) = tanh(mean(atanh(max(min(ra, 0.999), -0.999)), 1, 'omitnan'));
            gA = Ssel.group{find(strcmp(Ssel.animal, aList{a}), 1)};
            perAnimal(end+1) = struct('animal', aList{a}, 'group', gA, 'trialType', tt, ...
                'learnBin', lb, 'r', A(a,:), 'nSessions', sum(strcmp(Ssel.animal, aList{a}))); %#ok<AGROW>
        end

        mR = nan(1, numel(lags)); pR = nan(1, numel(lags));
        for li = 1:numel(lags)
            [mR(li), pR(li)] = signFlipMean(A(:, li));
        end
        [~, iPk] = max(abs(mR));
        stats = [stats; table(string(tt), string(lb), numel(aList), lagSec(iPk), mR(iPk), pR(iPk), ...
            {mR}, {pR}, ...
            'VariableNames', {'trialType','learnBin','nAnimals','peakLagSec','peakR','peakP','meanR','pByLag'})]; %#ok<AGROW>
    end
end
perAnimal = struct2table(perAnimal);

nA = max(stats.nAnimals);
fprintf('\n%s | window [%g %g] s | bins %g s | lags %+.2f to %+.2f s\n', ...
    sigLabel, opt.window(1), opt.window(2), dtBin, lagSec(1), lagSec(end));
fprintf('positive lag = LICKING LEADS DA | sign-flip p floor with %d animals: %.4f\n', nA, 2/2^nA);
disp(stats(:, {'trialType','learnBin','nAnimals','peakLagSec','peakR','peakP'}));

%% ---- figure ----
nRow = numel(trialTypes); nCol = max(1, numel(groupsToRun));
fig = figure('Color','w', 'Position', [100 100 380*nCol 320*nRow]);
cols = [0.80 0.80 0.91; 0.42 0.42 0.96; 0.05 0.05 0.86];
for t = 1:nRow
    for gI = 1:nCol
        lb = groupsToRun{gI};
        ax = subplot(nRow, nCol, (t-1)*nCol + gI); hold(ax, 'on');
        sel = strcmp(perAnimal.trialType, trialTypes{t});
        if useBin, sel = sel & perAnimal.learnBin == lb; end
        Pa = perAnimal(sel, :);
        if isempty(Pa), axis(ax,'off'); continue; end
        A = rowsOf(Pa.r, true(height(Pa), 1));
        for a = 1:size(A,1)
            c = [0 0.45 0.74]; if strcmp(Pa.group{a}, 'fast'), c = [0.85 0.33 0.10]; end
            plot(ax, lagSec, A(a,:), '-', 'Color', [c 0.35], 'LineWidth', 1);
        end
        mR = tanh(mean(atanh(max(min(A,0.999),-0.999)), 1, 'omitnan'));
        cc = [0 0 0]; if useBin, cc = cols(gI,:); end
        plot(ax, lagSec, mR, '-', 'Color', cc, 'LineWidth', 2.5);

        % mark lags significant across animals
        sRow = stats(strcmp(stats.trialType, trialTypes{t}) & stats.learnBin == string(lb), :);
        if ~isempty(sRow)
            pv = sRow.pByLag{1};
            sig = lagSec(pv < 0.05);
            if ~isempty(sig)
                yl = ylim(ax);
                plot(ax, sig, (yl(1) + 0.05*diff(yl))*ones(size(sig)), 's', ...
                    'MarkerSize', 4, 'MarkerFaceColor', cc, 'MarkerEdgeColor', 'none');
            end
        end
        xline(ax, 0, '--', 'Color', [0.5 0.5 0.5]);
        yline(ax, 0, '-',  'Color', [0.8 0.8 0.8]);
        xlabel(ax, 'lag (s)   \leftarrow DA leads    licking leads \rightarrow');
        ylabel(ax, sprintf('r (across trials), %s', upperType(trialTypes{t})));
        ttl = upperType(trialTypes{t});
        if useBin, ttl = sprintf('%s | %s', ttl, lb); end
        title(ax, ttl, 'FontWeight', 'normal');
        set(ax, 'TickDir','out'); grid(ax,'on'); box(ax,'off'); hold(ax,'off');
    end
end
sgtitle(fig, sprintf('%s vs lick count, trial-by-trial at each lag | window [%g %g] s', ...
    sigLabel, opt.window(1), opt.window(2)), 'FontWeight','bold','FontSize',11);

out = struct('perSession', perSession, 'perAnimal', perAnimal, 'stats', stats, ...
    'lagSec', lagSec, 'fig', fig);

if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d,'dir') ~= 7, mkdir(d); end
    f = fullfile(d, sprintf('DALickLagCorr_%s_%g-%g.pdf', sigTag, opt.window(1), opt.window(2)));
    set(fig,'InvertHardcopy','off');
    print(fig, f, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved:\n  %s\n', f);
end
end

%% ========================================================================
function M = rowsOf(col, mask)
% struct2table stores a same-length row-vector field as a NUMERIC matrix,
% but as a CELL column if any entry differs in length (e.g. a session that
% produced fewer lags). Accept either so the caller does not have to care.
if iscell(col)
    M = cell2mat(col(mask));
else
    M = col(mask, :);
end
end

%% ========================================================================
function [meanV, p] = signFlipMean(v)
v = v(isfinite(v));
n = numel(v);
if n < 2, meanV = NaN; p = NaN; return; end
z = atanh(max(min(v, 0.999), -0.999));
obs = mean(z);
signs = 1 - 2 * (dec2bin(0:2^n-1, n) - '0');
null = signs * z(:) / n;
p = sum(abs(null) >= abs(obs) - 1e-12) / size(signs,1);
meanV = tanh(obs);
end

function s = upperType(tt)
switch lower(tt)
    case 'cr', s = 'CR';
    case 'fa', s = 'FA';
    otherwise, s = [upper(tt(1)) tt(2:end)];
end
end
