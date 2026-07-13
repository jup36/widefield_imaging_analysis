function rez = computePerMouseLOOProjectionAndStats_correctTrials(anchorLOO, headerC, glmRezC, trIdC, expertHeaders_perMouse, varargin)
%computePerMouseLOOProjectionAndStats_correctTrials
% Project all sessions onto per-mouse LOO anchors, compute:
%   (1) Collapsed session metrics per time (per axis): muGo, muNoGo, diff, energy
%   (2) Permutation stats over session order (shared permutations across metrics)
%   (3) Trial-level omnibus KW across sessions (fast permutation-based)
%
% THIS VARIANT: muGo / muNoGo (and therefore diff, and optionally energy)
% are computed from CORRECT TRIALS ONLY by default:
%   muGo   <- mean projection over HIT trials   (not all Go trials)
%   muNoGo <- mean projection over CR trials    (not all NoGo trials)
%
% WHY
%   The original computePerMouseLOOProjectionAndStats pools ALL Go trials
%   (Hit + Miss) into muGo and ALL NoGo trials (CR + FA) into muNoGo, for
%   every session and every axis, uniformly. Because outcome mixture
%   changes systematically over learning (many FA/Miss early, mostly
%   CR/Hit once expert), any session-order trend (slope/Spearman) computed
%   on muGo/muNoGo/diff is confounded with the behavioral-outcome
%   trajectory, not just the neural one. Restricting to Hit/CR trials
%   holds the outcome category fixed across sessions, so a trend in the
%   correct-trial-only trajectory reflects the representation on
%   correctly-performed trials specifically.
%
% NEW BEHAVIOR (relative to the original function)
%   - goMask   is built from trId.hitI  instead of trId.goI   (by default)
%   - nogoMask is built from trId.crI   instead of trId.nogoI (by default)
%   - Legacy behavior is still available via 'TrialPolicy',"GoNoGoByGroup"
%     for A/B comparison against the original function.
%   - Each session's Hit-trial count and CR-trial count are recorded
%     (rezM.sessions.nGoUsed / nNoGoUsed). Sessions with fewer than
%     'MinTrialsPerCond' Hit trials OR fewer than 'MinTrialsPerCond' CR
%     trials are EXCLUDED from the trend/permutation stats (folded into
%     statSessMask, same mechanism as day4MarkC filtering) — NOT silently
%     included with a noisy few-trial mean. The collapsed muGo/muNoGo/diff
%     arrays still contain a value (possibly based on very few trials) for
%     every session for plotting purposes; rezM.sessions.sufficientTrials
%     tells you which sessions were trustworthy enough to be used for the
%     stats.
%   - 'EnergyTrialSubset' controls whether the "energy" collapsed metric
%     (mean squared projection, a non-directional metric) is computed over
%     ALL trials ("all", legacy behavior) or over the union of Hit and CR
%     trials only ("correctOnly", default) — since Miss/FA trials could
%     otherwise bias the overall energy trend independent of Go/NoGo
%     representation strength on correctly performed trials.
%
% REQUIRED
%   anchorLOO : output of buildPerMouseLOOAnchors_fromExperts
%               (or buildPerMouseLOOAnchors_fromExperts_correctTrials)
%   headerC, glmRezC, trIdC : (J x S) cells
%               trIdC{j,s} must have fields hitI/missI/crI/faI (and
%               goI/nogoI if TrialPolicy="GoNoGoByGroup" is used instead)
%   expertHeaders_perMouse  : (nMouse x nExpertMax) cells
%
% NAME-VALUE (common; identical to original unless noted)
%   WhichAxes        : "A" (default) | "Araw_ord"
%   DoCollapsed      : true (default)
%   DoOmnibusTrial   : true (default)
%
%   StatTypes        : {'slope','spearman'} (default)
%   nPerm            : 1000 (default)
%   FWER             : "both" (default)  % "max"|"cluster"|"both"|"none"
%   AlphaCluster     : 0.05 (default)
%   ClusterMassMode  : "sum" (default)   % "sum"|"sumabs"
%
%   day4MarkC        : [] (default) or Nx2 cell/string array
%                      e.g. {"m1045","121124"; ...}
%                      If non-empty, slope/spearman stats use only sessions
%                      with session date >= day4 date for that mouse.
%
%   NEW:
%   TrialPolicy        : "CorrectOnlyByGroup" (default) | "GoNoGoByGroup"
%                        Controls which trId fields define goMask/nogoMask.
%                        "CorrectOnlyByGroup" -> hitI / crI
%                        "GoNoGoByGroup"      -> goI  / nogoI (legacy)
%   MinTrialsPerCond   : 5 (default)
%                        Minimum Hit-trial count AND minimum CR-trial count
%                        (or Go/NoGo counts under legacy policy) required
%                        for a session to be included in trend/permutation
%                        stats. Sessions below threshold are excluded via
%                        statSessMask (AND'd with any day4MarkC filtering);
%                        their collapsed muGo/muNoGo/diff values are still
%                        computed and stored (for plotting) but flagged via
%                        rezM.sessions.sufficientTrials = false.
%   EnergyTrialSubset  : "correctOnly" (default) | "all"
%                        Which trials populate the "energy" collapsed
%                        metric. "correctOnly" = Hit union CR trials
%                        (or Go union NoGo under legacy policy); "all" =
%                        every trial regardless of outcome (legacy
%                        behavior).
%
%   SaveEachMouse    : true (default)
%   SaveDir          : 'Z:\Rodent Data\dualImaging_parkj\collectData\glmTDR_perMouseRez' (default)
%
%   Verbose          : true (default)
%   VerboseEvery     : 1 (default)
%   VerbosePerm      : true (default)
%   PermProgressPct  : 10 (default)

% -------------------- required checks --------------------
assert(iscell(headerC) && iscell(glmRezC) && iscell(trIdC), ...
    'headerC/glmRezC/trIdC must be cell.');
assert(isequal(size(headerC), size(glmRezC), size(trIdC)), ...
    'headerC/glmRezC/trIdC must be same size.');
assert(isstruct(anchorLOO) && isfield(anchorLOO,'perMouse'), ...
    'anchorLOO looks invalid.');
assert(iscell(expertHeaders_perMouse) || isstring(expertHeaders_perMouse), ...
    'expertHeaders_perMouse must be cell or string.');

% -------------------- defaults --------------------
opt = struct();
opt.WhichAxes       = "A";
opt.DoCollapsed     = true;
opt.DoOmnibusTrial  = false;

opt.StatTypes       = {'slope','spearman'};
opt.nPerm           = 1000;
opt.FWER            = "both";     % max|cluster|both|none
opt.AlphaCluster    = 0.05;
opt.ClusterMassMode = "sum";      % sum|sumabs

opt.day4MarkC       = [];

% NEW
opt.TrialPolicy       = "CorrectOnlyByGroup";  % "CorrectOnlyByGroup" | "GoNoGoByGroup"
opt.MinTrialsPerCond  = 5;
opt.EnergyTrialSubset = "correctOnly";         % "correctOnly" | "all"

opt.SaveEachMouse   = true;
opt.SaveDir         = 'Z:\Rodent Data\dualImaging_parkj\collectData\glmTDR_perMouseRez';

opt.Verbose         = true;
opt.VerboseEvery    = 1;
opt.VerbosePerm     = true;
opt.PermProgressPct = 10;

% -------------------- parse name-value (lightweight) --------------------
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0, 'Name-value args must come in pairs.');
    for i = 1:2:numel(varargin)
        k = char(string(varargin{i}));
        v = varargin{i+1};
        assert(isfield(opt, k), 'Unknown name-value: %s', k);
        opt.(k) = v;
    end
end

opt.WhichAxes       = string(opt.WhichAxes);
opt.FWER            = lower(string(opt.FWER));
opt.ClusterMassMode = lower(string(opt.ClusterMassMode));
opt.TrialPolicy       = string(opt.TrialPolicy);
opt.EnergyTrialSubset = lower(string(opt.EnergyTrialSubset));

if ischar(opt.StatTypes) || isstring(opt.StatTypes)
    opt.StatTypes = cellstr(string(opt.StatTypes));
end
opt.StatTypes = cellfun(@(s) lower(char(string(s))), opt.StatTypes, 'UniformOutput', false);

% determine trId field names for go-like / nogo-like conditions
switch opt.TrialPolicy
    case "CorrectOnlyByGroup"
        goFieldName   = 'hitI';
        nogoFieldName = 'crI';
    case "GoNoGoByGroup"
        goFieldName   = 'goI';
        nogoFieldName = 'nogoI';
    otherwise
        error('Unknown TrialPolicy: %s (must be "CorrectOnlyByGroup" or "GoNoGoByGroup")', opt.TrialPolicy);
end

% normalize day4 map
day4Map = [];
if ~isempty(opt.day4MarkC)
    day4Map = normalizeDay4MarkC_local_(opt.day4MarkC);
end

% create save dir if needed
if opt.SaveEachMouse
    saveDir = char(opt.SaveDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
end

% -------------------- trim trailing empties --------------------
[headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_local_(headerC, glmRezC, trIdC);

% -------------------- build header lookup --------------------
[uniqHdr, glmFirst, trFirst] = buildHeaderLookup_local_(headerC, glmRezC, trIdC);

% -------------------- normalize expert header matrix --------------------
hdrMat = normalizeHeaderMatrix_local_(expertHeaders_perMouse);
nMouse = size(hdrMat,1);

% -------------------- output shell --------------------
rez = struct();
rez.opt = opt;
rez.perMouse = cell(nMouse, 1);

% ==================== main loop over mice ====================
for m = 1:nMouse

    % -------- mouse id --------
    mouseId = "";
    if isfield(anchorLOO.perMouse,'mouseIds') && numel(anchorLOO.perMouse.mouseIds) >= m ...
            && ~isempty(anchorLOO.perMouse.mouseIds{m})
        mouseId = string(anchorLOO.perMouse.mouseIds{m});
    end
    if strlength(mouseId)==0
        mouseId = inferMouseIdFromRow_local_(hdrMat(m,:));
    end
    if strlength(mouseId)==0
        mouseId = "mouse" + string(m);
    end

    if opt.Verbose && (mod(m,opt.VerboseEvery)==0 || m==1 || m==nMouse)
        fprintf('[Proj+Stats:correctTrials] mouse %d/%d (%s)\n', m, nMouse, mouseId);
    end

    % -------- collect sessions for this mouse --------
    [sessHdr, sessGlm, sessTr] = collectMouseSessions_local_(mouseId, uniqHdr, glmFirst, trFirst);
    nSess = numel(sessHdr);

    rezM = struct();
    rezM.mouseId = char(mouseId);
    rezM.sessions = struct();
    rezM.sessions.headers = cellstr(sessHdr(:));
    rezM.sessions.nSess = nSess;

    if nSess == 0
        warning('Mouse %s: no sessions found in headerC lookup. Skipping.', char(mouseId));
        rez.perMouse{m} = rezM;
        continue;
    end

    % -------- session date parsing (used by day4 filter too) --------
    [sessDt, sessDtOk] = parseHeaderDatesOnly_local_(sessHdr);
    rezM.sessions.date = sessDt;

    % -------- expert membership --------
    expertRow = string(hdrMat(m,:));
    expertRow = expertRow(strlength(expertRow)>0);
    isExpertSess = ismember(sessHdr, expertRow);
    rezM.sessions.isExpert = isExpertSess;

    % -------- axesFull / axesLOO / kMap --------
    axesFull = [];
    if isfield(anchorLOO.perMouse,'axesFull') && numel(anchorLOO.perMouse.axesFull) >= m
        axesFull = anchorLOO.perMouse.axesFull{m};
    end
    axesLOO = [];
    if isfield(anchorLOO.perMouse,'axesLOO')
        axesLOO = anchorLOO.perMouse.axesLOO(m,:);
    end
    kMap = [];
    if isfield(anchorLOO.perMouse,'kMapHeaderToIdx') && numel(anchorLOO.perMouse.kMapHeaderToIdx) >= m
        kMap = anchorLOO.perMouse.kMapHeaderToIdx{m};
    end

    [Afull, axisNames] = getAxisMatrix_local_(axesFull, opt.WhichAxes);
    if isempty(Afull)
        warning('Mouse %s: axesFull empty. Skipping.', char(mouseId));
        rez.perMouse{m} = rezM;
        continue;
    end
    nAxis = size(Afull,1);

    rezM.axes = struct();
    rezM.axes.which = char(opt.WhichAxes);
    rezM.axes.nAxis = nAxis;
    rezM.axes.names = axisNames;

    % -------- time base --------
    timeVec = sessGlm{1}.decBins.time(:);
    nTime = numel(timeVec);
    rezM.time = timeVec;

    % -------- allocate collapsed metrics --------
    if opt.DoCollapsed
        muGo   = nan(nSess, nTime, nAxis, 'single');
        muNoGo = nan(nSess, nTime, nAxis, 'single');
        diffMN = nan(nSess, nTime, nAxis, 'single');
        energy = nan(nSess, nTime, nAxis, 'single');
    end

    % -------- NEW: per-session trial-count bookkeeping --------
    nGoUsed   = nan(nSess,1);
    nNoGoUsed = nan(nSess,1);

    % -------- project each session --------
    for s = 1:nSess
        gr = sessGlm{s};
        tr = sessTr{s};

        % choose axes for this session
        Ause = Afull;
        if isExpertSess(s)
            if ~isempty(kMap) && isa(kMap,'containers.Map') && isKey(kMap, char(sessHdr(s)))
                kk = kMap(char(sessHdr(s)));
                if ~isempty(axesLOO) && kk >= 1 && kk <= numel(axesLOO) && ~isempty(axesLOO{kk})
                    [Aloo, ~] = getAxisMatrix_local_(axesLOO{kk}, opt.WhichAxes);
                    if ~isempty(Aloo)
                        Ause = Aloo;
                    end
                end
            end
        end

        [Z, N, nW] = projectSession_local_(gr, Ause); % [trial x time x axis]
        assert(nW == nTime, 'Time bins mismatch in session %s', char(sessHdr(s)));

        % ---- CORRECT-TRIALS-ONLY (by default) go/nogo masks ----
        goMask   = local_makeMask_from_trId_(tr, goFieldName,   N);
        nogoMask = local_makeMask_from_trId_(tr, nogoFieldName, N);

        nGo   = sum(goMask);
        nNoGo = sum(nogoMask);
        nGoUsed(s)   = nGo;
        nNoGoUsed(s) = nNoGo;

        if opt.DoCollapsed
            if nGo >= 1
                Zgo = Z(goMask,:,:);
                muGo(s,:,:) = single(squeeze(mean(Zgo, 1, 'omitnan')));
            end
            if nNoGo >= 1
                Znogo = Z(nogoMask,:,:);
                muNoGo(s,:,:) = single(squeeze(mean(Znogo, 1, 'omitnan')));
            end
            diffMN(s,:,:) = muGo(s,:,:) - muNoGo(s,:,:); % NaN if either side missing

            switch opt.EnergyTrialSubset
                case "correctonly"
                    energyMask = goMask | nogoMask;
                    if any(energyMask)
                        energy(s,:,:) = single(squeeze(mean(Z(energyMask,:,:).^2, 1, 'omitnan')));
                    end
                case "all"
                    energy(s,:,:) = single(squeeze(mean(Z.^2, 1, 'omitnan')));
                otherwise
                    error('Unknown EnergyTrialSubset: %s', opt.EnergyTrialSubset);
            end
        end

        if opt.Verbose && (nGo < opt.MinTrialsPerCond || nNoGo < opt.MinTrialsPerCond)
            fprintf('  [lowTrials] %s | nGo(%s)=%d nNoGo(%s)=%d (Min=%d) -> excluded from trend stats\n', ...
                char(sessHdr(s)), goFieldName, nGo, nogoFieldName, nNoGo, opt.MinTrialsPerCond);
        end
    end

    if opt.DoCollapsed
        rezM.collapsed = struct();
        rezM.collapsed.muGo   = muGo;
        rezM.collapsed.muNoGo = muNoGo;
        rezM.collapsed.diff   = diffMN;
        rezM.collapsed.energy = energy;
    end

    rezM.sessions.nGoUsed   = nGoUsed;
    rezM.sessions.nNoGoUsed = nNoGoUsed;
    sufficientTrials = (nGoUsed >= opt.MinTrialsPerCond) & (nNoGoUsed >= opt.MinTrialsPerCond);
    rezM.sessions.sufficientTrials = sufficientTrials;

    % -------- determine sessions used for trend stats --------
    statSessMask = true(nSess,1);
    day4Date = NaT;

    if ~isempty(day4Map)
        hit = find(strcmpi(day4Map.mouseId, mouseId), 1, 'first');
        if ~isempty(hit)
            day4Date = day4Map.date(hit);
            if ~sessDtOk
                error('Mouse %s: could not parse session dates, but day4MarkC filtering was requested.', char(mouseId));
            end
            statSessMask = sessDt >= day4Date;
        else
            error('Mouse %s not found in day4MarkC.', char(mouseId));
        end
    end

    % NEW: AND in the trial-sufficiency gate
    statSessMask = statSessMask & sufficientTrials;

    rezM.sessions.statSessMask = statSessMask;
    rezM.sessions.day4Date     = day4Date;

    if opt.Verbose
        nExcludedByTrials = sum(~sufficientTrials);
        fprintf('  [trialGate] %s | %d/%d sessions excluded from stats for insufficient %s/%s trials (Min=%d)\n', ...
            char(mouseId), nExcludedByTrials, nSess, goFieldName, nogoFieldName, opt.MinTrialsPerCond);
        if ~isempty(day4Map)
            fprintf('  [day4] %s | keep %d/%d sessions from %s onward (before trial gate)\n', ...
                char(mouseId), sum(sessDt >= day4Date), nSess, char(string(day4Date, 'MMddyy')));
        end
        fprintf('  [stats-sessions] %s | final nSessStat = %d/%d (day4 AND trialGate)\n', ...
            char(mouseId), sum(statSessMask), nSess);
    end

    % -------- permutation stats over sessions --------
    rezM.stats = struct();

    if opt.DoCollapsed && opt.nPerm > 0
        Y = struct();
        Y.muGo   = muGo;
        Y.muNoGo = muNoGo;
        Y.diff   = diffMN;
        Y.energy = energy;

        metricNames = fieldnames(Y);
        reshapeStat = @(v) reshape(v, [nTime, nAxis])'; % 1 x (nTime*nAxis) -> [nAxis x nTime]

        % apply day4 + trial-gate filtering for stats
        statIdx = find(statSessMask);
        nSessStat = numel(statIdx);
        rezM.sessions.nSessStat = nSessStat;
        rezM.sessions.statIdx   = statIdx(:);

        if nSessStat < 2
            warning('Mouse %s: fewer than 2 sessions remain after stat-session filtering (day4 + trial gate). Skipping trend stats.', char(mouseId));
        else
            x = (1:nSessStat)';

            for stI = 1:numel(opt.StatTypes)
                statType = string(opt.StatTypes{stI});

                if opt.Verbose
                    fprintf('  [stats] %s | nPerm=%d | FWER=%s | share-perm across metrics | nSessStat=%d\n', ...
                        char(statType), opt.nPerm, char(opt.FWER), nSessStat);
                end

                % Vectorize filtered metrics to [nSessStat x (nTime*nAxis)]
                Y2 = struct();
                for mi = 1:numel(metricNames)
                    mn = metricNames{mi};
                    Yf = Y.(mn)(statIdx,:,:); % filtered sessions only
                    Y2.(mn) = double(reshape(Yf, [nSessStat, nTime*nAxis]));
                end

                obsAll = struct();

                if statType == "slope"
                    xc = x - mean(x);
                    denomX = max(sum(xc.^2), eps);
                    sdX = std(x);

                    obsAll_SD = struct();

                    for mi = 1:numel(metricNames)
                        mn = metricNames{mi};

                        Ymat = Y2.(mn);                    % [nSessStat x cols]
                        vObs = (xc' * Ymat) / denomX;     % raw slope

                        sdY = std(Ymat, 0, 1, 'omitnan');
                        vObsSD = vObs .* (sdX ./ max(sdY, eps));

                        obsAll.(mn)    = reshapeStat(vObs);
                        obsAll_SD.(mn) = reshapeStat(vObsSD);
                    end

                elseif statType == "spearman"
                    rx = tiedrank(x);
                    rx = rx - mean(rx);
                    denomX_obs = sqrt(max(sum(rx.^2), eps));

                    spearman_rYc    = struct();
                    spearman_denomY = struct();

                    for mi = 1:numel(metricNames)
                        mn = metricNames{mi};

                        Ymat = Y2.(mn);
                        rY = zeros(size(Ymat), 'double');
                        for c = 1:size(Ymat,2)
                            rY(:,c) = tiedrank(Ymat(:,c));
                        end
                        rYc = rY - mean(rY,1);
                        denomY = sqrt(max(sum(rYc.^2,1), eps));

                        spearman_rYc.(mn)    = rYc;
                        spearman_denomY.(mn) = denomY;

                        vObs = (rx' * rYc) ./ max(denomX_obs .* denomY, eps);
                        obsAll.(mn) = reshapeStat(vObs);
                    end
                else
                    error('Unknown StatType: %s', char(statType));
                end

                % Allocate permutation stats
                permAll = struct();
                for mi = 1:numel(metricNames)
                    mn = metricNames{mi};
                    permAll.(mn) = zeros(opt.nPerm, nAxis, nTime, 'single');
                end

                reportEvery = max(1, round(opt.nPerm * (opt.PermProgressPct/100)));
                if opt.VerbosePerm
                    fprintf('    [perm] report every %d perms (~%d%%)\n', reportEvery, opt.PermProgressPct);
                end

                for p = 1:opt.nPerm
                    idx = randperm(nSessStat);
                    xP  = x(idx);

                    if statType == "slope"
                        xcP = xP - mean(xP);
                        denomXP = max(sum(xcP.^2), eps);

                        for mi = 1:numel(metricNames)
                            mn = metricNames{mi};
                            Ymat = Y2.(mn);
                            vP = (xcP' * Ymat) / denomXP;   % raw slope
                            permAll.(mn)(p,:,:) = single(reshapeStat(vP));
                        end

                    elseif statType == "spearman"
                        rxP = tiedrank(xP);
                        rxP = rxP - mean(rxP);
                        denomXP = sqrt(max(sum(rxP.^2), eps));

                        for mi = 1:numel(metricNames)
                            mn = metricNames{mi};
                            rYc    = spearman_rYc.(mn);
                            denomY = spearman_denomY.(mn);

                            vP = (rxP' * rYc) ./ max(denomXP .* denomY, eps);
                            permAll.(mn)(p,:,:) = single(reshapeStat(vP));
                        end
                    end

                    if opt.VerbosePerm && (p==1 || p==opt.nPerm || mod(p,reportEvery)==0)
                        pct = floor(100*p/opt.nPerm);
                        fprintf('      perm %d/%d (%d%%)\n', p, opt.nPerm, pct);
                    end
                end

                % Compute p-values + store
                for mi = 1:numel(metricNames)
                    mn = metricNames{mi};

                    statObs  = double(obsAll.(mn));   % [nAxis x nTime]
                    statPerm = double(permAll.(mn));  % [nPerm x nAxis x nTime]

                    out = struct();
                    out.stat = single(statObs);
                    out.nSessStat = nSessStat;
                    out.statIdx   = statIdx(:);
                    out.day4Date  = day4Date;

                    if statType == "slope"
                        out.statSD = single(double(obsAll_SD.(mn)));
                    end

                    % ---------- two-sided uncorrected ----------
                    p_unc = ones(nAxis, nTime, 'single');
                    for a = 1:nAxis
                        permAT = squeeze(statPerm(:,a,:)); % [nPerm x nTime]
                        obsAT  = statObs(a,:);             % [1 x nTime]
                        p_unc(a,:) = single((1 + sum(abs(permAT) >= abs(obsAT), 1)) ./ (opt.nPerm + 1));
                    end
                    out.p_unc = p_unc;

                    % ---------- one-sided uncorrected (direction matched to observed sign) ----------
                    p1_unc = ones(nAxis, nTime, 'single');
                    for a = 1:nAxis
                        permAT = squeeze(statPerm(:,a,:));   % [nPerm x nTime]
                        obsAT  = statObs(a,:);               % [1 x nTime]

                        ptmp = ones(1, nTime);
                        posI = obsAT > 0;
                        negI = obsAT < 0;
                        zerI = obsAT == 0;

                        if any(posI)
                            ptmp(posI) = (1 + sum(permAT(:,posI) >= obsAT(posI), 1)) ./ (opt.nPerm + 1);
                        end
                        if any(negI)
                            ptmp(negI) = (1 + sum(permAT(:,negI) <= obsAT(negI), 1)) ./ (opt.nPerm + 1);
                        end
                        if any(zerI)
                            ptmp(zerI) = 1;
                        end

                        p1_unc(a,:) = single(ptmp);
                    end
                    out.p1_unc = p1_unc;

                    % ---------- two-sided max-FWER ----------
                    if opt.FWER == "max" || opt.FWER == "both"
                        p_max = ones(nAxis, nTime, 'single');
                        for a = 1:nAxis
                            permAT  = squeeze(statPerm(:,a,:));
                            maxNull = max(abs(permAT), [], 2);
                            obsAT   = abs(statObs(a,:));
                            p_max(a,:) = single((1 + sum(maxNull >= obsAT, 1)) ./ (opt.nPerm + 1));
                        end
                        out.p_max = p_max;
                    end

                    % ---------- one-sided max-FWER (direction matched to observed sign) ----------
                    if opt.FWER == "max" || opt.FWER == "both"
                        p1_max = ones(nAxis, nTime, 'single');
                        for a = 1:nAxis
                            permAT   = squeeze(statPerm(:,a,:));   % [nPerm x nTime]
                            maxNullP = max(permAT, [], 2);         % upper-tail max null
                            minNullP = min(permAT, [], 2);         % lower-tail min null
                            obsAT    = statObs(a,:);               % [1 x nTime]

                            ptmp = ones(1, nTime);
                            posI = obsAT > 0;
                            negI = obsAT < 0;
                            zerI = obsAT == 0;

                            if any(posI)
                                ptmp(posI) = (1 + sum(maxNullP >= obsAT(posI), 1)) ./ (opt.nPerm + 1);
                            end
                            if any(negI)
                                ptmp(negI) = (1 + sum(minNullP <= obsAT(negI), 1)) ./ (opt.nPerm + 1);
                            end
                            if any(zerI)
                                ptmp(zerI) = 1;
                            end

                            p1_max(a,:) = single(ptmp);
                        end
                        out.p1_max = p1_max;
                    end

                    % ---------- two-sided cluster-FWER ----------
                    if opt.FWER == "cluster" || opt.FWER == "both"
                        p_cluster   = ones(nAxis, nTime, 'single');
                        thr_cluster = nan(nAxis, nTime, 'single');
                        clusterInfo = cell(nAxis,1);

                        for a = 1:nAxis
                            permAT = squeeze(statPerm(:,a,:));
                            obsAT  = statObs(a,:).';
                            thr    = prctile(abs(permAT), 100*(1-opt.AlphaCluster), 1).';

                            [pC, info] = clusterFWER_local_(obsAT, permAT, thr, opt.ClusterMassMode);

                            p_cluster(a,:)   = single(pC(:)).';
                            thr_cluster(a,:) = single(thr(:)).';
                            clusterInfo{a}   = info;
                        end

                        out.p_cluster   = p_cluster;
                        out.thr_cluster = thr_cluster;
                        out.cluster     = clusterInfo;
                    end

                    % ---------- one-sided cluster-FWER (direction matched to observed sign) ----------
                    if opt.FWER == "cluster" || opt.FWER == "both"
                        p1_cluster   = ones(nAxis, nTime, 'single');
                        thr1_cluster = nan(nAxis, nTime, 'single');
                        clusterInfo1 = cell(nAxis,1);

                        for a = 1:nAxis
                            permAT = squeeze(statPerm(:,a,:));   % [nPerm x nTime]
                            obsAT  = statObs(a,:).';             % [nTime x 1]

                            posI = obsAT > 0;
                            negI = obsAT < 0;

                            pAxis   = ones(nTime,1);
                            thrAxis = nan(nTime,1);
                            infoAxis = struct();

                            % positive-direction clusters
                            if any(posI)
                                thrPos = prctile(permAT, 100*(1-opt.AlphaCluster), 1).';
                                [pPos, infoPos] = clusterFWER_local_positive_(obsAT, permAT, thrPos, opt.ClusterMassMode);

                                thrAxis(posI) = thrPos(posI);
                                pAxis(posI)   = pPos(posI);
                                infoAxis.pos  = infoPos;
                            else
                                infoAxis.pos = [];
                            end

                            % negative-direction clusters
                            if any(negI)
                                thrNeg = prctile(permAT, 100*(opt.AlphaCluster), 1).';
                                [pNeg, infoNeg] = clusterFWER_local_negative_(obsAT, permAT, thrNeg, opt.ClusterMassMode);

                                thrAxis(negI) = thrNeg(negI);
                                pAxis(negI)   = pNeg(negI);
                                infoAxis.neg  = infoNeg;
                            else
                                infoAxis.neg = [];
                            end

                            % zeros stay p=1
                            p1_cluster(a,:)   = single(pAxis(:)).';
                            thr1_cluster(a,:) = single(thrAxis(:)).';
                            clusterInfo1{a}   = infoAxis;
                        end

                        out.p1_cluster   = p1_cluster;
                        out.thr1_cluster = thr1_cluster;
                        out.cluster1     = clusterInfo1;
                    end

                    rezM.stats.(char(statType)).(mn) = out;
                end
            end
        end
    end

    % -------- write back --------
    rez.perMouse{m} = rezM;

    % -------- save per mouse --------
    if opt.SaveEachMouse
        saveDir = char(opt.SaveDir);
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end

        mmddyy = char(datetime('today','Format','MMddyy'));
        fn = sprintf('%s_glmTDR_perMouseRez_correctTrials_%s.mat', char(mouseId), mmddyy);
        fpath = fullfile(saveDir, fn);

        rezMouse = rezM; %#ok<NASGU>
        optUsed  = opt;  %#ok<NASGU>
        try
            save(fpath, 'rezMouse', 'optUsed', '-v7.3');
            if opt.Verbose
                fprintf('  [save] %s\n', fpath);
            end
        catch ME
            warning('Save failed for %s: %s', fpath, ME.message);
        end
    end

end % for m

end % function


%% ======================================================================
% Helpers (local subfunctions) — unchanged from original unless noted
% ======================================================================

function day4 = normalizeDay4MarkC_local_(day4In)
if isstring(day4In)
    day4In = cellstr(day4In);
end
assert(iscell(day4In), 'day4MarkC must be cell or string.');
assert(size(day4In,2) == 2, 'day4MarkC must be Nx2: {mouseId, "MMDDYY"; ...}.');

n = size(day4In,1);
day4 = struct();
day4.mouseId = strings(n,1);
day4.date    = NaT(n,1);

for i = 1:n
    mid = string(day4In{i,1});
    dstr = string(day4In{i,2});
    assert(strlength(mid)>0, 'day4MarkC row %d has empty mouseId.', i);
    assert(strlength(dstr)>0, 'day4MarkC row %d has empty date.', i);

    try
        d = datetime(char(dstr), 'InputFormat','MMddyy');
    catch
        error('day4MarkC row %d has invalid date string: %s', i, char(dstr));
    end

    day4.mouseId(i) = mid;
    day4.date(i)    = d;
end
end

function [sessHdr, sessGlm, sessTr] = collectMouseSessions_local_(mouseId, uniqHdr, glmFirst, trFirst)
mouseId = string(mouseId);
mask = contains(uniqHdr, mouseId);
sessHdr = uniqHdr(mask);
sessGlm = glmFirst(mask);
sessTr  = trFirst(mask);

[dtKey, ok] = parseHeaderDatetimeWithSuffix_local_(sessHdr);
if ok
    [~, ord] = sort(dtKey, 'ascend');
    sessHdr = sessHdr(ord);
    sessGlm = sessGlm(ord);
    sessTr  = sessTr(ord);
end
end

function [dtOnly, ok] = parseHeaderDatesOnly_local_(hdrS)
hdrS = string(hdrS(:));
n = numel(hdrS);
dtOnly = NaT(n,1);
ok = true;
for i = 1:n
    h = char(hdrS(i));
    tok = regexp(h, '_(\d{6})(?:-(\d+))?$', 'tokens', 'once');
    if isempty(tok)
        ok = false;
        return;
    end
    mmddyy = tok{1};
    try
        dtOnly(i) = datetime(mmddyy, 'InputFormat','MMddyy');
    catch
        ok = false;
        return;
    end
end
end

function [Z, N, nW] = projectSession_local_(glmRez, Ause)
Yz = glmRez.Yz;                          % [M x K]
timeVec = glmRez.decBins.time(:);
nW = numel(timeVec);
M = size(Yz,1);
assert(rem(M,nW)==0, 'M not divisible by nW');
N = M / nW;

Zrows = Yz * Ause';                      % [M x nAxis]
nAxis = size(Ause,1);
Z = reshape(Zrows, [N, nW, nAxis]);     % [trial x time x axis]
end

function m = local_makeMask_from_trId_(trId, field, N)
m = true(N,1);
if isempty(trId) || ~isstruct(trId) || ~isfield(trId, field) || isempty(trId.(field))
    warning('local_makeMask_from_trId_:MissingField', ...
        'trId.%s missing/empty; using all trials for this condition.', field);
    return;
end
x = trId.(field);
if islogical(x)
    x = x(:);
    if numel(x)==N
        m = x;
    else
        warning('local_makeMask_from_trId_:BadLength', ...
            'trId.%s length %d != N=%d; using all trials for this condition.', field, numel(x), N);
    end
else
    idx = unique(round(x(:)));
    idx = idx(idx>=1 & idx<=N);
    m = false(N,1);
    m(idx) = true;
end
end

function out = omnibusKW_perm_fast_overTime_local_(trialCell, nTime, nAxis, opt)
% Fast permutation-based Kruskal-Wallis omnibus test across sessions.
% (Unchanged from original; DoOmnibusTrial defaults to false.)

nSess = numel(trialCell);

out = struct();
out.stat = nan(nTime, nAxis, 'single');
out.df   = nan(nTime, nAxis, 'single');

if nSess < 2
    out.p_unc = nan(nTime, nAxis, 'single');
    return;
end

ranksCell = cell(nAxis, nTime);
gCell     = cell(nAxis, nTime);
NCell     = nan(nAxis, nTime);
dfCell    = nan(nAxis, nTime);

for a = 1:nAxis
    for t = 1:nTime
        xAll = [];
        gAll = [];

        for s = 1:nSess
            Xs = trialCell{s};
            if isempty(Xs)
                continue;
            end

            v = double(Xs(:,t,a));
            v = v(isfinite(v));
            if isempty(v)
                continue;
            end

            xAll = [xAll; v]; %#ok<AGROW>
            gAll = [gAll; s*ones(numel(v),1)]; %#ok<AGROW>
        end

        if numel(xAll) < 2 || numel(unique(gAll)) < 2
            continue;
        end

        ranksCell{a,t} = tiedrank(xAll);
        gCell{a,t}     = gAll;
        NCell(a,t)     = numel(xAll);
        dfCell(a,t)    = numel(unique(gAll)) - 1;

        out.df(t,a)   = single(dfCell(a,t));
        out.stat(t,a) = single(kw_chi2_fast_local_(ranksCell{a,t}, gAll, NCell(a,t)));
    end
end

permStat = nan(opt.nPerm, nAxis, nTime, 'single');

reportEvery = max(1, round(opt.nPerm * (opt.PermProgressPct/100)));
if opt.VerbosePerm
    fprintf('    [omnibus perm] report every %d perms (~%d%%)\n', reportEvery, opt.PermProgressPct);
end

for p = 1:opt.nPerm
    for a = 1:nAxis
        for t = 1:nTime
            ranks = ranksCell{a,t};
            gAll  = gCell{a,t};

            if isempty(ranks)
                continue;
            end

            gPerm = gAll(randperm(numel(gAll)));
            permStat(p,a,t) = single(kw_chi2_fast_local_(ranks, gPerm, NCell(a,t)));
        end
    end

    if opt.VerbosePerm && (p==1 || p==opt.nPerm || mod(p,reportEvery)==0)
        pct = floor(100*p/opt.nPerm);
        fprintf('      perm %d/%d (%d%%)\n', p, opt.nPerm, pct);
    end
end

out.p_unc = ones(nTime, nAxis, 'single');

for a = 1:nAxis
    obsAT  = double(out.stat(:,a))';
    permAT = squeeze(double(permStat(:,a,:)));

    validT = isfinite(obsAT);
    if any(validT)
        ptmp = ones(1, nTime);
        ptmp(validT) = (1 + sum(permAT(:,validT) >= obsAT(validT), 1)) ./ (opt.nPerm + 1);
        out.p_unc(:,a) = single(ptmp(:));
    end
end

if opt.FWER == "max" || opt.FWER == "both"
    out.p_max = ones(nTime, nAxis, 'single');

    for a = 1:nAxis
        obsAT  = double(out.stat(:,a))';
        permAT = squeeze(double(permStat(:,a,:)));
        maxNull = max(permAT, [], 2);

        validT = isfinite(obsAT);
        if any(validT)
            ptmp = ones(1, nTime);
            ptmp(validT) = (1 + sum(maxNull >= obsAT(validT), 1)) ./ (opt.nPerm + 1);
            out.p_max(:,a) = single(ptmp(:));
        end
    end
end

if opt.FWER == "cluster" || opt.FWER == "both"
    out.p_cluster   = ones(nTime, nAxis, 'single');
    out.thr_cluster = nan(nTime, nAxis, 'single');
    out.cluster     = cell(nAxis,1);

    for a = 1:nAxis
        obsAT  = double(out.stat(:,a));
        permAT = squeeze(double(permStat(:,a,:)));

        validT = isfinite(obsAT);
        if ~any(validT)
            continue;
        end

        thr = nan(nTime,1);
        thr(validT) = prctile(permAT(:,validT), 100*(1-opt.AlphaCluster), 1)';

        [pC, info] = clusterFWER_local_positive_(obsAT, permAT, thr, opt.ClusterMassMode);

        out.p_cluster(:,a)   = single(pC(:));
        out.thr_cluster(:,a) = single(thr(:));
        out.cluster{a}       = info;
    end
end
end

function [A, names] = getAxisMatrix_local_(axesStruct, whichAxes)
A = [];
names = {};
if isempty(axesStruct) || ~isstruct(axesStruct)
    return;
end
whichAxes = string(whichAxes);

if whichAxes == "A"
    if isfield(axesStruct,'A') && ~isempty(axesStruct.A)
        A = axesStruct.A;
    end
elseif whichAxes == "Araw_ord"
    if isfield(axesStruct,'Araw_ord') && ~isempty(axesStruct.Araw_ord)
        A = axesStruct.Araw_ord;
    end
else
    error('WhichAxes must be "A" or "Araw_ord"');
end

if isfield(axesStruct,'names') && ~isempty(axesStruct.names) && whichAxes=="A"
    names = axesStruct.names;
elseif isfield(axesStruct,'names_ord') && ~isempty(axesStruct.names_ord) && whichAxes=="Araw_ord"
    names = axesStruct.names_ord;
elseif isfield(axesStruct,'names_ord') && ~isempty(axesStruct.names_ord)
    names = axesStruct.names_ord;
else
    names = arrayfun(@(i) sprintf('axis_%d', i), 1:size(A,1), 'UniformOutput', false);
end
end

function hdrMat = normalizeHeaderMatrix_local_(hdrIn)
if isempty(hdrIn)
    hdrMat = cell(0,0);
    return;
end
if isstring(hdrIn)
    hdrIn = cellstr(hdrIn);
end
hdrMat = cell(size(hdrIn));
for i = 1:numel(hdrIn)
    x = hdrIn{i};
    if isempty(x)
        hdrMat{i} = [];
    else
        s = char(string(x));
        if strlength(string(s))==0
            hdrMat{i} = [];
        else
            hdrMat{i} = s;
        end
    end
end
end

function mouseId = inferMouseIdFromRow_local_(hdrRow)
mouseId = "";
for k = 1:numel(hdrRow)
    h = hdrRow{k};
    if isempty(h)
        continue;
    end
    tok = regexp(string(h), '(m\d{3,5})', 'tokens', 'once');
    if ~isempty(tok)
        mouseId = string(tok{1});
        return;
    end
end
end

function [uniqHdr, glmFlatFirst, trFlatFirst] = buildHeaderLookup_local_(headerC, glmRezC, trIdC)
hdrFlat = headerC(:);
glmFlat = glmRezC(:);
trFlat  = trIdC(:);

isHdr = ~cellfun(@isempty, hdrFlat);
hdrFlat = hdrFlat(isHdr);
glmFlat = glmFlat(isHdr);
trFlat  = trFlat(isHdr);

hdrFlatS = string(hdrFlat);
[uniqHdr, ~, ic] = unique(hdrFlatS, 'stable');

if numel(uniqHdr) < numel(hdrFlatS)
    counts = accumarray(ic, 1);
    dup = uniqHdr(counts > 1);
    warning('Duplicate headers in headerC; using first occurrence. Example: %s', string(dup(1)));
end

glmFlatFirst = cell(numel(uniqHdr),1);
trFlatFirst  = cell(numel(uniqHdr),1);
for i = 1:numel(uniqHdr)
    ii = find(hdrFlatS == uniqHdr(i), 1, 'first');
    glmFlatFirst{i} = glmFlat{ii};
    trFlatFirst{i}  = trFlat{ii};
end
end

function [dtKey, ok] = parseHeaderDatetimeWithSuffix_local_(hdrS)
hdrS = string(hdrS(:));
n = numel(hdrS);
dtKey = NaT(n,1);
ok = true;
for i = 1:n
    h = char(hdrS(i));
    tok = regexp(h, '_(\d{6})(?:-(\d+))?$', 'tokens', 'once');
    if isempty(tok)
        ok = false;
        return;
    end
    mmddyy = tok{1};
    suf = 0;
    if numel(tok) >= 2 && ~isempty(tok{2})
        suf = str2double(tok{2});
        if ~isfinite(suf)
            suf = 0;
        end
    end
    d0 = datetime(mmddyy, 'InputFormat','MMddyy');
    dtKey(i) = d0 + seconds(suf);
end
end

function [headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_local_(headerC, glmRezC, trIdC)
[J,~] = size(headerC);
Smax = max([size(headerC,2), size(glmRezC,2), size(trIdC,2)]);
if size(headerC,2) < Smax
    headerC(:,end+1:Smax) = {[]};
end
if size(glmRezC,2) < Smax
    glmRezC(:,end+1:Smax) = {[]};
end
if size(trIdC,2) < Smax
    trIdC(:,end+1:Smax) = {[]};
end

keepLast = 0;
for s = 1:Smax
    anyNonEmpty = any(~cellfun(@isempty, headerC(:,s))) || ...
                  any(~cellfun(@isempty, glmRezC(:,s)))  || ...
                  any(~cellfun(@isempty, trIdC(:,s)));
    if anyNonEmpty
        keepLast = s;
    end
end

if keepLast == 0
    headerC = cell(J,0);
    glmRezC = cell(J,0);
    trIdC   = cell(J,0);
else
    headerC = headerC(:,1:keepLast);
    glmRezC = glmRezC(:,1:keepLast);
    trIdC   = trIdC(:,1:keepLast);
end
end

function [p_clust, info] = clusterFWER_local_(statObs, statPerm, thr, massMode)
% Two-sided cluster-FWER for signed stats.

nPerm = size(statPerm,1);
nTime = numel(statObs);

sup = abs(statObs) > thr;
cl = findClusters_local_(sup);

obsMass = zeros(numel(cl),1);
for c = 1:numel(cl)
    idx = cl{c};
    switch massMode
        case "sum"
            obsMass(c) = sum(statObs(idx));
        case "sumabs"
            obsMass(c) = sum(abs(statObs(idx)));
        otherwise
            error('Unknown ClusterMassMode: %s', massMode);
    end
end
obsMassAbs = abs(obsMass);

maxNull = zeros(nPerm,1);
thrCol = thr(:);

for p = 1:nPerm
    st = statPerm(p,:)';
    supP = abs(st) > thrCol;
    clP = findClusters_local_(supP);
    if isempty(clP)
        maxNull(p) = 0;
        continue;
    end
    mP = zeros(numel(clP),1);
    for c = 1:numel(clP)
        idx = clP{c};
        switch massMode
            case "sum"
                mP(c) = sum(st(idx));
            case "sumabs"
                mP(c) = sum(abs(st(idx)));
        end
    end
    maxNull(p) = max(abs(mP));
end

p_clust = ones(nTime,1);
pClusterEach = nan(numel(cl),1);
for c = 1:numel(cl)
    pval = (1 + sum(maxNull >= obsMassAbs(c))) / (nPerm + 1);
    pClusterEach(c) = pval;
    p_clust(cl{c}) = pval;
end

info = struct();
info.nClusters = numel(cl);
info.clusters  = cl;
info.massObs   = obsMass;
info.pCluster  = pClusterEach;
end

function cl = findClusters_local_(mask)
mask = mask(:);
d = diff([false; mask; false]);
starts = find(d==1);
ends   = find(d==-1)-1;
cl = cell(numel(starts),1);
for i = 1:numel(starts)
    cl{i} = starts(i):ends(i);
end
end

function [p_clust, info] = clusterFWER_local_positive_(statObs, statPerm, thr, massMode)
% One-sided / positive-only cluster-FWER for positive-valued stats.

nPerm = size(statPerm,1);
nTime = numel(statObs);

sup = statObs > thr;
cl = findClusters_local_(sup);

obsMass = zeros(numel(cl),1);
for c = 1:numel(cl)
    idx = cl{c};
    switch massMode
        case "sum"
            obsMass(c) = sum(statObs(idx));
        case "sumabs"
            obsMass(c) = sum(statObs(idx)); % same for positive stats
        otherwise
            error('Unknown ClusterMassMode: %s', massMode);
    end
end

maxNull = zeros(nPerm,1);

for p = 1:nPerm
    st = statPerm(p,:)';
    supP = st > thr;
    clP = findClusters_local_(supP);

    if isempty(clP)
        maxNull(p) = 0;
        continue;
    end

    mP = zeros(numel(clP),1);
    for c = 1:numel(clP)
        idx = clP{c};
        switch massMode
            case "sum"
                mP(c) = sum(st(idx));
            case "sumabs"
                mP(c) = sum(st(idx));
        end
    end
    maxNull(p) = max(mP);
end

p_clust = ones(nTime,1);
pClusterEach = nan(numel(cl),1);

for c = 1:numel(cl)
    pval = (1 + sum(maxNull >= obsMass(c))) / (nPerm + 1);
    pClusterEach(c) = pval;
    p_clust(cl{c}) = pval;
end

info = struct();
info.nClusters = numel(cl);
info.clusters  = cl;
info.massObs   = obsMass;
info.pCluster  = pClusterEach;
end

function [p_clust, info] = clusterFWER_local_negative_(statObs, statPerm, thr, massMode)
% One-sided / negative-only cluster-FWER for signed stats.
% statObs:  [nTime x 1]
% statPerm: [nPerm x nTime]
% thr:      [nTime x 1] lower threshold (typically alpha quantile)

nPerm = size(statPerm,1);
nTime = numel(statObs);

sup = statObs < thr;
cl = findClusters_local_(sup);

obsMass = zeros(numel(cl),1);
for c = 1:numel(cl)
    idx = cl{c};
    switch massMode
        case "sum"
            obsMass(c) = sum(statObs(idx));          % negative mass
        case "sumabs"
            obsMass(c) = sum(abs(statObs(idx)));     % magnitude-only
        otherwise
            error('Unknown ClusterMassMode: %s', massMode);
    end
end

minNull = zeros(nPerm,1);

for p = 1:nPerm
    st = statPerm(p,:)';
    supP = st < thr;
    clP = findClusters_local_(supP);

    if isempty(clP)
        switch massMode
            case "sum"
                minNull(p) = 0;
            case "sumabs"
                minNull(p) = 0;
        end
        continue;
    end

    mP = zeros(numel(clP),1);
    for c = 1:numel(clP)
        idx = clP{c};
        switch massMode
            case "sum"
                mP(c) = sum(st(idx));              % negative values
            case "sumabs"
                mP(c) = sum(abs(st(idx)));
        end
    end

    switch massMode
        case "sum"
            minNull(p) = min(mP);                 % most negative cluster
        case "sumabs"
            minNull(p) = max(mP);                 % largest magnitude
    end
end

p_clust = ones(nTime,1);
pClusterEach = nan(numel(cl),1);

for c = 1:numel(cl)
    switch massMode
        case "sum"
            pval = (1 + sum(minNull <= obsMass(c))) / (nPerm + 1);
        case "sumabs"
            pval = (1 + sum(minNull >= obsMass(c))) / (nPerm + 1);
    end
    pClusterEach(c) = pval;
    p_clust(cl{c}) = pval;
end

info = struct();
info.nClusters = numel(cl);
info.clusters  = cl;
info.massObs   = obsMass;
info.pCluster  = pClusterEach;
end

function chi2 = kw_chi2_fast_local_(ranks, g, N)
% Fast Kruskal-Wallis chi-square statistic from precomputed ranks.

if isempty(ranks) || isempty(g) || N < 2
    chi2 = NaN;
    return;
end

groups = unique(g);
k = numel(groups);
if k < 2
    chi2 = NaN;
    return;
end

sumTerm = 0;
for i = 1:k
    gi = groups(i);
    idx = (g == gi);
    ni = sum(idx);
    if ni == 0
        continue;
    end
    Ri = sum(ranks(idx));
    sumTerm = sumTerm + (Ri.^2) / ni;
end

chi2 = (12 / (N * (N + 1))) * sumTerm - 3 * (N + 1);
end