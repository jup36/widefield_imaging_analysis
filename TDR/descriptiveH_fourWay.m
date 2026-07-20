function HstatC = descriptiveH_fourWay(HsY3d, Htime, trI, varargin)
%DESCRIPTIVEH_FOURWAY  Compute descriptive statistics (mean & SEM) of motif
% activity (H), split into the four trial-OUTCOME types: Hit, Miss, CR, FA.
%
% This is the four-outcome extension of descriptiveH.m (which only splits
% by stimulus IDENTITY: Go vs NoGo, i.e. Hit+Miss pooled vs CR+FA pooled).
% Splitting by outcome instead of stimulus identity lets you see, motif by
% motif, how the raw temporal activation trace differs between correct and
% incorrect trials of the same stimulus type -- e.g. whether a
% NoGo-related motif's activity on FA trials looks more like its Hit-trial
% profile than its CR-trial profile. This is intended as an intuitive,
% single-motif companion to the GLM-TDR projection-score analyses: H is
% the raw (z-scored) motif activation itself, not a GLM-fitted or
% axis-projected quantity.
%
% HstatC = descriptiveH_fourWay(HsY3d, Htime, trI, 'Name', value, ...)
%
% INPUTS
%   HsY3d : [N x K x T] motif activity array (trial x motif x time),
%           typically Hs.Y3 from stack_trials_H(tbytDat_hAligned, 'zscore', true)
%   Htime : [1 x T] time vector (e.g. Hs.winCtrs). Stored for reference
%           only; not used internally.
%   trI   : struct with logical masks hitI, missI, crI, faI (each [N x 1]
%           or [1 x N]), e.g. from trialTypeInfoAuditoryGngTbytDat(tbytDat)
%
% NAME-VALUE
%   'MinTrialsPerCond' : 3 (default). Minimum trial count required to
%                        compute mean/SEM for a given motif x condition.
%                        Below this, mean/SEM are set to NaN(1,T) rather
%                        than silently computed from a handful of trials
%                        -- consistent with the "no silent fallback"
%                        principle used elsewhere in this pipeline.
%                        nTrial is still recorded regardless, so you can
%                        see why a condition was excluded.
%   'Verbose'          : false (default)
%
% OUTPUT
%   HstatC : 1 x K cell array; HstatC{k} has fields:
%       .motifIdx
%       .mean.hit / .mean.miss / .mean.cr / .mean.fa   -- each [1 x T] (or NaN(1,T))
%       .sem.hit  / .sem.miss  / .sem.cr  / .sem.fa
%       .nTrial.hit / .nTrial.miss / .nTrial.cr / .nTrial.fa
%
% DEPENDENCIES
%   meanstdsem.m
%       Expected signature: [m, sd, se] = meanstdsem(data)
%       Operates along dim 1 (rows = trials).
%
% EXAMPLE
%   Hs = stack_trials_H(tbytDat_hAligned, 'zscore', true);
%   trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);
%   HstatC = descriptiveH_fourWay(Hs.Y3, Hs.winCtrs, trI, 'MinTrialsPerCond', 3);
%
%   % Plot motif #5, all four trial types
%   figure; hold on;
%   plot(Hs.winCtrs, HstatC{5}.mean.hit,  'b', 'LineWidth', 2);
%   plot(Hs.winCtrs, HstatC{5}.mean.miss, 'b--');
%   plot(Hs.winCtrs, HstatC{5}.mean.cr,   'r', 'LineWidth', 2);
%   plot(Hs.winCtrs, HstatC{5}.mean.fa,   'r--');
%   legend('Hit','Miss','CR','FA');
%
% See also: descriptiveH (Go/NoGo two-way version), stack_trials_H,
%           trialTypeInfoAuditoryGngTbytDat

ip = inputParser;
ip.FunctionName = mfilename;
ip.addParameter('MinTrialsPerCond', 3, @(x) isnumeric(x) && isscalar(x) && x>=0);
ip.addParameter('Verbose', false, @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});
opt = ip.Results;

assert(ndims(HsY3d) == 3, 'HsY3d must be [N x K x T] (trial x motif x time).'); %#ok<ISMAT>
[N, K, T] = size(HsY3d);

condFields     = {'hit','miss','cr','fa'};
condMaskFields = struct('hit','hitI','miss','missI','cr','crI','fa','faI');

% -------------------- validate + coerce trial masks --------------------
maskC = struct();
for c = 1:numel(condFields)
    cf = condFields{c};
    mf = condMaskFields.(cf);
    assert(isfield(trI, mf), 'trI is missing required field "%s".', mf);

    m = trI.(mf);
    m = logical(m(:));
    assert(numel(m) == N, ...
        'trI.%s length (%d) does not match N (%d, from HsY3d dim 1).', mf, numel(m), N);

    maskC.(cf) = m;
end

% -------------------- per-motif descriptive stats --------------------
HstatC = cell(1, K);

for k = 1:K
    Hstat = struct();
    Hstat.motifIdx = k;
    Hstat.mean   = struct();
    Hstat.sem    = struct();
    Hstat.nTrial = struct();

    Xk = squeeze(HsY3d(:, k, :));   % [N x T]
    if isvector(Xk)
        % Guard against squeeze collapsing a singleton T or N dimension
        Xk = reshape(Xk, N, T);
    end

    for c = 1:numel(condFields)
        cf = condFields{c};
        m  = maskC.(cf);
        nC = sum(m);

        Hstat.nTrial.(cf) = nC;

        if nC >= opt.MinTrialsPerCond
            [mC, ~, seC] = meanstdsem(Xk(m, :));
            Hstat.mean.(cf) = mC;
            Hstat.sem.(cf)  = seC;
        else
            if opt.Verbose
                fprintf(['[descriptiveH_fourWay] motif %d, cond "%s": ' ...
                    'nTrial=%d < MinTrialsPerCond=%d -> NaN\n'], ...
                    k, cf, nC, opt.MinTrialsPerCond);
            end
            Hstat.mean.(cf) = nan(1, T);
            Hstat.sem.(cf)  = nan(1, T);
        end
    end

    HstatC{1, k} = Hstat;
end

end
