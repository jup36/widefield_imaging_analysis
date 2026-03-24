function out = tdrSingleSessionDiag(glmRez, trI, varargin)
% ============================= SYNOPSIS =============================
% tdrSingleSessionDiag.m
%
% PURPOSE
%   Run GLM-based Targeted Dimensionality Reduction (GLM-TDR) on a *single*
%   session (one glmRez struct) and perform a quick diagnostic check that
%   Go vs NoGo trials separate along a user-specified axis (e.g., GoToneOn_1).
%   Optionally, visualize a raw observed signal (one column of Yz or Ybig)
%   in the same Go/NoGo image layout to sanity-check trial indexing and
%   stacking conventions.
%
% WHAT IT DOES (HIGH-LEVEL)
%   1) Recompute per-session GLM-TDR axes from glmRez using targetedDimRed_from_glmRez
%        - Uses OrthMode="none" so the axis set is *not* Gram–Schmidt orthonormalized.
%        - Axis selection is done by name match in tdr.names_ord (pre-GS ordered axes).
%   2) Select Go vs NoGo trial sets from trI
%        - Default: Go vs NoGo uses trI.goI and trI.nogoI
%        - If CorrectTrialsOnly=true: Hit vs CR uses trI.hitI and trI.crI
%        - Trial selectors can be numeric indices or logical masks; validated against N.
%   3) Project the session’s activity matrix (Yz or Ybig) onto the specified axis
%        - Axis vector is taken from tdr.Araw_ord(axisIdx,:) and re-normalized.
%        - Projection is done in row-space: zRows = Yproj * axis'
%   4) Reshape the projected vector into [N x nW] using the TIME-MAJOR convention
%        - TIME-MAJOR means rows are stacked as:
%            [timeBin1 trials 1..N; timeBin2 trials 1..N; ...; timeBin(nW) trials 1..N]
%        - nW is inferred from numel(glmRez.decBins.time)
%        - N is inferred as M / nW where M = size(Yproj,1)
%   5) Quantify a simple separation metric
%        - score_timeMajor = mean over time of |mean(Go) - mean(NoGo)|
%        - This is not a classifier; it’s just a quick sanity metric.
%
% OPTIONAL: “OBSERVED” VIEW (raw signal sanity check)
%   If ShowObserved=true (default), the function also:
%     - Extracts one feature column from ObservedWhichY (Yz or Ybig),
%       e.g., ObservedFeatureIdx=1 for motif 1
%     - Reshapes it into the same [N x nW] time-major trial-by-time matrix
%     - Plots Go vs NoGo images and reports an analogous separation score
%   This is useful for detecting broken trial indexing or inconsistent stacking,
%   because it bypasses axis construction and uses a raw signal.
%
% VISUALIZATION OUTPUTS
%   If MakeFigures=true, generates:
%     - Figure 1: imagesc of projected axis activity
%         left = Go (or Hit), right = NoGo (or CR)
%     - Figure 2 (optional): imagesc of observed raw feature activity
%         left = Go/Hit, right = NoGo/CR
%
% IMPORTANT NEW FEATURE: CLIMIT (consistent color scaling)
%   Name-value: 'climit', [] (default) or [cmin cmax]
%     - If provided, applies caxis(climit) to ALL imagesc panels (projection + observed).
%     - This prevents misleading comparisons caused by per-panel autoscaling.
%     - Strongly recommended when visually comparing Go vs NoGo.
%
% INPUTS / EXPECTATIONS
%   glmRez must contain:
%     - glmRez.decBins.time (defines nW)
%     - glmRez.Yz (and/or glmRez.Ybig) with size [M x K]
%     - fields required by targetedDimRed_from_glmRez (beta, X_design, muX, sdX, group, ...)
%   trI must contain either:
%     - goI / nogoI  (default comparison) OR
%     - hitI / crI   (CorrectTrialsOnly=true)
%
% OUTPUT (out struct)
%   out.tdr               : full per-session TDR output (from targetedDimRed_from_glmRez)
%   out.axisName, axisIdx : which axis was used
%   out.N, out.nW, out.K  : inferred dimensions
%   out.goI, out.nogoI    : resolved numeric trial indices used
%   out.Z_timeMajor       : [N x nW] projected activity for the axis
%   out.score_timeMajor   : simple separation score
%   out.obs (optional)    : observed raw feature matrices + score (if ShowObserved=true)
%
% ====================================================================

% -------------------- parse --------------------
p = inputParser;
p.addRequired('glmRez', @(s)isstruct(s) && ~isempty(s));
p.addRequired('trI',    @(s)isstruct(s) && ~isempty(s));

p.addParameter('ProjectWhichY', "Yz", @(s)ischar(s)||isstring(s));
p.addParameter('AxisName', "GoToneOn_1", @(s)ischar(s)||isstring(s));
p.addParameter('CorrectTrialsOnly', false, @(x)islogical(x)&&isscalar(x));

p.addParameter('ShowObserved', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ObservedWhichY', "Yz", @(s)ischar(s)||isstring(s));
p.addParameter('ObservedFeatureIdx', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));

p.addParameter('SignFix', "maxabs", @(s)ischar(s)||isstring(s));
p.addParameter('Eps', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('MakeFigures', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

% NEW
p.addParameter('climit', [], @(x)isnumeric(x) && (isempty(x) || numel(x)==2));

p.parse(glmRez, trI, varargin{:});
opt = p.Results;

% -------------------- TDR --------------------
tdr = targetedDimRed_from_glmRez(glmRez, ...
    'OrthMode', "none", ...
    'SignFix',  opt.SignFix, ...
    'Eps',      opt.Eps, ...
    'ProjectWhichY', char(string(opt.ProjectWhichY)), ...
    'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
    'MultiDimK', opt.MultiDimK, ...
    'PriorityNames', cellstr(string(opt.PriorityNames)));

% -------------------- trial selection --------------------
if opt.CorrectTrialsOnly
    goSel   = trI.hitI;
    nogoSel = trI.crI;
    labGo   = 'Hit';
    labNoGo = 'CR';
else
    goSel   = trI.goI;
    nogoSel = trI.nogoI;
    labGo   = 'Go';
    labNoGo = 'NoGo';
end

nW = numel(glmRez.decBins.time);
[Yproj, N, K] = getY_(glmRez, opt.ProjectWhichY, nW);

goI   = normalizeTrialSelector_(goSel, N, 'goI/hitI');
nogoI = normalizeTrialSelector_(nogoSel, N, 'nogoI/crI');

% -------------------- axis --------------------
axisName = char(string(opt.AxisName));
namesOrd = cellstr(string(tdr.names_ord));
idxAxis = find(strcmpi(namesOrd, axisName), 1);

Arow = tdr.Araw_ord(idxAxis,:);
Arow = Arow ./ max(norm(Arow), opt.Eps);

zRows = Yproj * Arow';
Z_timeMajor = reshape(zRows, [N, nW]);

score_time = mean(abs(mean(Z_timeMajor(goI,:),1) - mean(Z_timeMajor(nogoI,:),1)));

% -------------------- observed --------------------
obs = [];
if opt.ShowObserved
    [Yobs, Nobs, Kobs] = getY_(glmRez, opt.ObservedWhichY, nW);
    fIdx = opt.ObservedFeatureIdx;

    yRows = Yobs(:, fIdx);
    Y_timeMajor = reshape(yRows, [N, nW]);

    obs_score_time = mean(abs(mean(Y_timeMajor(goI,:),1) - mean(Y_timeMajor(nogoI,:),1)));

    obs = struct();
    obs.whichY = char(string(opt.ObservedWhichY));
    obs.featureIdx = fIdx;
    obs.Y_timeMajor = Y_timeMajor;
    obs.score_timeMajor = obs_score_time;
end

% -------------------- plotting --------------------
if opt.MakeFigures

    % -------- Projection --------
    figure('Name', sprintf('Projection | %s', axisName));
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    imagesc(Z_timeMajor(goI,:));
    title(sprintf('%s (n=%d)', labGo, numel(goI)));
    xlabel('time'); ylabel('trials'); colorbar;
    if ~isempty(opt.climit), caxis(opt.climit); end

    nexttile;
    imagesc(Z_timeMajor(nogoI,:));
    title(sprintf('%s (n=%d)', labNoGo, numel(nogoI)));
    xlabel('time'); ylabel('trials'); colorbar;
    if ~isempty(opt.climit), caxis(opt.climit); end

    % -------- Observed --------
    if opt.ShowObserved
        figure('Name', sprintf('Observed %s(:,%d)', obs.whichY, obs.featureIdx));
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

        nexttile;
        imagesc(obs.Y_timeMajor(goI,:));
        title(sprintf('%s (n=%d)', labGo, numel(goI)));
        xlabel('time'); ylabel('trials'); colorbar;
        if ~isempty(opt.climit), caxis(opt.climit); end

        nexttile;
        imagesc(obs.Y_timeMajor(nogoI,:));
        title(sprintf('%s (n=%d)', labNoGo, numel(nogoI)));
        xlabel('time'); ylabel('trials'); colorbar;
        if ~isempty(opt.climit), caxis(opt.climit); end
    end
end

% -------------------- output --------------------
out = struct();
out.tdr = tdr;
out.axisName = axisName;
out.axisIdx  = idxAxis;
out.N = N;
out.nW = nW;
out.K = K;
out.goI = goI;
out.nogoI = nogoI;
out.Z_timeMajor = Z_timeMajor;
out.score_timeMajor = score_time;
out.obs = obs;

end