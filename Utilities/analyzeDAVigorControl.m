function out = analyzeDAVigorControl(rezC, sessT, varargin)
%ANALYZEDAVIGORCONTROL
%   Does the growth of DA across sessions survive controlling for how much
%   the animal licks on each trial? Two-stage: fit within each animal, then
%   test the resulting coefficients across animals.
%
%   STAGE 1, per animal -- pool that animal's trials across its sessions and
%   fit two ordinary least-squares models to the per-trial DA scalar:
%       (a) DA ~ sessOrder
%       (b) DA ~ sessOrder + lick
%   Keep b_uncontrolled and b_controlled, the sessOrder slope from each.
%   Pooling across sessions is what makes the question answerable: sessOrder
%   barely varies WITHIN a session, so a per-session fit would have almost no
%   leverage on it.
%
%   STAGE 2, across animals -- exact sign-flip permutation on the per-animal
%   slopes (all 2^n assignments; 512 for 9 mice). Three tests:
%       b_uncontrolled  : is there a session-order effect at all?
%       b_controlled    : does it survive the lick control?
%       b_controlled - b_uncontrolled : did controlling change it?
%   The animal is the unit of inference, the test is exact, and its
%   resolution is visible (p floor = 2/2^n).
%
%   HOW TO READ THE RESULT -- both directions are weaker than they look
%   ------------------------------------------------------------------
%   SURVIVES: rules out the simplest vigor account, that DA tracks how much
%     the animal licks in the window. It does NOT rule out vigor on
%     dimensions lick count misses (force, tongue kinematics, posture). And
%     lick count is a noisy proxy for vigor: measurement error in a covariate
%     biases its own coefficient toward zero and leaves variance for
%     sessOrder to absorb, so "survives" overstates its case.
%   DOES NOT SURVIVE: this is NOT evidence that DA is merely movement. Lick
%     vigor is plausibly downstream of the same learning driving DA -- a
%     mediator, not a nuisance -- and conditioning on a mediator removes the
%     effect of interest by construction.
%   For that reason the model-free companions matter more than the
%   regression: 'matched' (early vs late trials at comparable lick counts)
%   and the Hit/FA contrast (vigorous licking with and without reward).
%
%   out = analyzeDAVigorControl(rezC, sessT, ...)
%
% NAME-VALUE
%   'signal'       : 'global' (default) or motif number 1..K
%   'trialTypes'   : {'hit','cr'} -- analysed separately
%   'window'       : 'resp' (default) | 'cue'  -- which per-trial scalar
%   'minTrials'    : 30  -- an animal needs this many usable trials
%   'zscoreWithin' : true -- z-score DA, lick and sessOrder within animal
%                    before fitting, so slopes are comparable across animals
%                    of different scale (they become standardized betas)
%   'matched'      : true -- also run the matched-lick comparison
%   'matchTol'     : 1   -- licks; early/late trials pair within this
%   'fastLearners' : {'m1044','m1045','m1092','m1094'}
%   'slowLearners' : {'m1048','m1049','m1613','m1859','m1873'}
%   'figSaveDir'   : ''
%
% OUTPUT
%   .perAnimal  one row per animal x trialType: b_uncontrolled, b_controlled,
%               b_lick, nTrials, nSessions, matchedDiff
%   .stats      one row per trialType: mean slopes and sign-flip p's
%   .fig

%% ---- parse ----
p = inputParser;
p.addParameter('signal', 'global', @(x) (ischar(x)||isstring(x)) || (isnumeric(x)&&isscalar(x)));
p.addParameter('trialTypes', {'hit','cr'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('window', 'resp', @(s) any(strcmpi(string(s), ["resp","cue"])));
p.addParameter('minTrials', 30, @(x) isnumeric(x) && isscalar(x) && x >= 10);
p.addParameter('zscoreWithin', true, @(x) islogical(x) && isscalar(x));
p.addParameter('matched', true, @(x) islogical(x) && isscalar(x));
p.addParameter('matchTol', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('fastLearners', {'m1044','m1045','m1092','m1094'}, @iscell);
p.addParameter('slowLearners', {'m1048','m1049','m1613','m1859','m1873'}, @iscell);
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

trialTypes = cellstr(opt.trialTypes);
win = lower(char(string(opt.window)));
switch win
    case 'resp', daF = 'daResp'; lkF = 'lickResp';
    case 'cue',  daF = 'daCue';  lkF = 'lickCue';
end

if isnumeric(opt.signal)
    sigRow = opt.signal + 1;                     % row 1 is global
    sigLabel = sprintf('motif %d DA', opt.signal);
    sigTag = sprintf('motif%02d', opt.signal);
else
    assert(strcmpi(opt.signal,'global'), 'signal must be ''global'' or a motif number.');
    sigRow = 1; sigLabel = 'global DA'; sigTag = 'global';
end

%% ---- gather per-trial data, pooled within animal ----
animals = unique(string(sessT.animal), 'stable');
rows = struct('animal',{},'group',{},'trialType',{},'b_uncontrolled',{},'b_controlled',{}, ...
    'b_lick',{},'nTrials',{},'nSessions',{},'matchedDiff',{},'matchedN',{});

for a = 1:numel(animals)
    aId = char(animals(a));
    if ismember(aId, opt.fastLearners), g = 'fast';
    elseif ismember(aId, opt.slowLearners), g = 'slow';
    else, continue;
    end
    ia = find(string(sessT.animal) == animals(a));

    for t = 1:numel(trialTypes)
        tt = trialTypes{t};
        DA = []; LK = []; SO = []; nSess = 0;
        for i = ia(:)'
            r = rezC{sessT.a(i), sessT.s(i)};
            if ~isfield(r, 'trial')
                continue;   % session predates the per-trial arrays
            end
            m = r.trial.type == tt;
            if ~any(m), continue; end
            nSess = nSess + 1;
            DA = [DA, double(r.trial.(daF)(sigRow, m))];        %#ok<AGROW>
            LK = [LK, double(r.trial.(lkF)(m))];                %#ok<AGROW>
            SO = [SO, repmat(sessT.sessOrder(i), 1, sum(m))];   %#ok<AGROW>
        end

        ok = isfinite(DA) & isfinite(LK) & isfinite(SO);
        DA = DA(ok)'; LK = LK(ok)'; SO = SO(ok)';
        if numel(DA) < opt.minTrials || numel(unique(SO)) < 3 || std(LK) == 0
            continue;
        end

        if opt.zscoreWithin
            % Standardize within animal so the per-animal slopes are on one
            % scale before being averaged across animals.
            DAf = (DA - mean(DA)) / std(DA);
            LKf = (LK - mean(LK)) / std(LK);
            SOf = (SO - mean(SO)) / std(SO);
        else
            DAf = DA; LKf = LK; SOf = SO;
        end

        b0 = [ones(numel(SOf),1), SOf] \ DAf;                   % DA ~ sessOrder
        b1 = [ones(numel(SOf),1), SOf, LKf] \ DAf;              % DA ~ sessOrder + lick

        % matched-lick comparison: early vs late halves by session order,
        % restricted to trials whose lick counts pair within matchTol.
        mDiff = NaN; mN = 0;
        if opt.matched
            [mDiff, mN] = matchedLickDiff(DA, LK, SO, opt.matchTol);
        end

        rows(end+1) = struct('animal', aId, 'group', g, 'trialType', tt, ...
            'b_uncontrolled', b0(2), 'b_controlled', b1(2), 'b_lick', b1(3), ...
            'nTrials', numel(DA), 'nSessions', nSess, ...
            'matchedDiff', mDiff, 'matchedN', mN); %#ok<AGROW>
    end
end
assert(~isempty(rows), 'No animal had enough usable trials -- has the batch been rerun with r.trial?');
perAnimal = struct2table(rows);

%% ---- stage 2: sign-flip across animals ----
stats = table();
for t = 1:numel(trialTypes)
    tt = trialTypes{t};
    R = perAnimal(strcmp(perAnimal.trialType, tt), :);
    [m0, p0] = signFlipMean(R.b_uncontrolled);
    [m1, p1] = signFlipMean(R.b_controlled);
    [mD, pD] = signFlipMean(R.b_controlled - R.b_uncontrolled);
    [mL, pL] = signFlipMean(R.b_lick);
    [mM, pM] = signFlipMean(R.matchedDiff);
    stats = [stats; table(string(tt), height(R), m0, p0, m1, p1, mD, pD, mL, pL, mM, pM, ...
        'VariableNames', {'trialType','nAnimals', ...
            'mean_b_uncontrolled','p_uncontrolled', ...
            'mean_b_controlled','p_controlled', ...
            'mean_b_change','p_change', ...
            'mean_b_lick','p_lick', ...
            'mean_matchedDiff','p_matched'})]; %#ok<AGROW>
end

nA = max(stats.nAnimals);
fprintf('\n%s | %s window | two-stage (pooled within animal, sign-flip across)\n', sigLabel, win);
fprintf('sign-flip p floor with %d animals: %.4f\n', nA, 2/2^nA);
disp(stats);
fprintf(['Read p_controlled as the answer to "does the session-order effect survive?"; ' ...
         'p_matched is the model-free version.\n']);

%% ---- figure ----
fig = figure('Color','w', 'Position', [100 100 420*numel(trialTypes) 380]);
for t = 1:numel(trialTypes)
    tt = trialTypes{t};
    R = perAnimal(strcmp(perAnimal.trialType, tt), :);
    ax = subplot(1, numel(trialTypes), t); hold(ax, 'on');
    isF = strcmp(R.group, 'fast');
    for i = 1:height(R)
        c = [0 0.45 0.74]; if isF(i), c = [0.85 0.33 0.10]; end
        plot(ax, [1 2], [R.b_uncontrolled(i), R.b_controlled(i)], '-', ...
            'Color', [c 0.5], 'LineWidth', 1.1);
        scatter(ax, [1 2], [R.b_uncontrolled(i), R.b_controlled(i)], 45, c, 'filled');
    end
    errorbar(ax, [1 2], [mean(R.b_uncontrolled,'omitnan'), mean(R.b_controlled,'omitnan')], ...
        [semL(R.b_uncontrolled), semL(R.b_controlled)], 'ko-', ...
        'MarkerFaceColor','k','LineWidth',1.6,'CapSize',10);
    yline(ax, 0, '--', 'Color', [0.5 0.5 0.5]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'DA~session', 'DA~session+lick'}, 'XLim', [0.7 2.3]);
    ylabel(ax, 'session-order slope');
    s = stats(strcmp(stats.trialType, tt), :);
    title(ax, sprintf('%s | p = %.3g -> %.3g (n = %d mice)', upper(tt), ...
        s.p_uncontrolled, s.p_controlled, s.nAnimals), 'FontWeight','normal');
    set(ax, 'TickDir','out'); grid(ax,'on'); box(ax,'off'); hold(ax,'off');
end
sgtitle(fig, sprintf('%s, %s window | session-order slope before and after controlling for lick count', ...
    sigLabel, win), 'FontWeight','bold','FontSize',11);

out = struct('perAnimal', perAnimal, 'stats', stats, 'fig', fig);

if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d,'dir') ~= 7, mkdir(d); end
    f = fullfile(d, sprintf('DAVigorControl_%s_%s.pdf', sigTag, win));
    set(fig,'InvertHardcopy','off');
    print(fig, f, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved:\n  %s\n', f);
end
end

%% ========================================================================
function [d, nPair] = matchedLickDiff(DA, LK, SO, tol)
% Model-free companion: split trials at the animal's median session order,
% then compare late vs early DA only among trials with COMPARABLE lick
% counts. Greedy nearest-lick pairing without replacement. Because it never
% regresses DA on lick, it avoids both the mediator problem and any
% assumption that the DA-lick relation is linear.
d = NaN; nPair = 0;
mid = median(SO);
iE = find(SO <= mid); iL = find(SO > mid);
if isempty(iE) || isempty(iL), return; end

usedL = false(numel(iL),1);
dd = [];
for q = 1:numel(iE)
    cand = find(~usedL & abs(LK(iL) - LK(iE(q))) <= tol);
    if isempty(cand), continue; end
    [~, k] = min(abs(LK(iL(cand)) - LK(iE(q))));
    j = cand(k);
    usedL(j) = true;
    dd(end+1) = DA(iL(j)) - DA(iE(q)); %#ok<AGROW>
end
if isempty(dd), return; end
d = mean(dd); nPair = numel(dd);
end

%% ========================================================================
function [meanV, p] = signFlipMean(v)
% Exact sign-flip permutation over all 2^n assignments.
v = v(isfinite(v));
n = numel(v);
if n < 2, meanV = NaN; p = NaN; return; end
obs = mean(v);
signs = 1 - 2 * (dec2bin(0:2^n-1, n) - '0');
null = signs * v(:) / n;
p = sum(abs(null) >= abs(obs) - 1e-12) / size(signs,1);
meanV = obs;
end

function s = semL(v)
v = v(isfinite(v));
if numel(v) < 2, s = NaN; else, s = std(v)/sqrt(numel(v)); end
end