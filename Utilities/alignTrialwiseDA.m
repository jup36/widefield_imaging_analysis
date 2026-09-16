%% ========================================================================
function alignedC = alignTrialwiseDA(DAglobalC, tbytDat, tint)
% Interpolate each trial's raw DA trace onto tint (seconds from evtOn).
% No extrapolation: bins outside the trial's own coverage stay NaN, so a
% short trial is visibly short rather than padded with edge values.
nTrials = numel(tbytDat);
nTime = numel(tint);
alignedC = cell(2, nTrials);
for tr = 1:nTrials
    alignedC{1, tr} = NaN(1, nTime);
    alignedC{2, tr} = tint;
end

for tr = 1:nTrials
    y = DAglobalC{1, tr};
    tAbs = DAglobalC{2, tr};
    if isempty(y) || isempty(tAbs), continue; end

    evt = getNumericScalar(tbytDat(tr).evtOn);
    if ~isfinite(evt), continue; end

    tRel = tAbs - evt;
    ok = isfinite(tRel) & isfinite(y);
    if sum(ok) < 2, continue; end
    tRel = tRel(ok); y = y(ok);

    [tRel, si] = sort(tRel, 'ascend');  y = y(si);
    [tRel, ui] = unique(tRel, 'stable'); y = y(ui);
    if numel(tRel) < 2, continue; end

    inRange = tint >= tRel(1) & tint <= tRel(end);
    if ~any(inRange), continue; end

    alignedC{1, tr}(inRange) = interp1(tRel, y, tint(inRange), 'linear');
end
end
