function HestAll = motif_event_estimates( ...
        beta, X_names, muX, sdX, Xz, ...
        toneOn_Go_tbyb, toneOff_Go_tbyb, ...
        toneOn_NoGo_tbyb, toneOff_NoGo_tbyb, ...
        basisToneOn, basisToneOff, ...
        winCtrs, trI, varargin)
% MOTIF_EVENT_ESTIMATES  Event-locked estimated h (z-units) from β and bases.
%
% HestAll = motif_event_estimates( ...
%       beta, X_names, muX, sdX, Xz, ...
%       toneOn_Go_tbyb, toneOff_Go_tbyb, ...
%       toneOn_NoGo_tbyb, toneOff_NoGo_tbyb, ...
%       basisToneOn, basisToneOff, ...
%       winCtrs, trI, ...
%       'ReturnTrials', false)
%
% INPUTS (key ones)
%   beta     : [P x K] ridge coefficients (for standardized X and Yz)
%   X_names  : 1xP cellstr
%   muX      : [1 x P] column means used to z-score X
%   sdX      : [1 x P] column stds  used to z-score X
%   Xz       : [(N*nW) x P] full standardized design matrix (bin-major)
%
%   toneOn_*/toneOff_*_tbyb, basisToneOn/basisToneOff, winCtrs, trI :
%       as before (same ones used to build X).
%
% OUTPUT (per motif k)
%   HestAll{k}.mean.toneOnGo,  etc   : event-specific means (z-units)
%   HestAll{k}.mean.fullGo/NoGo      : full predicted H means (z-units)
%   HestAll{k}.sem.(...)             : SEMs for each mean
%   HestAll{k}.trial.(...)           : per-trial matrices if ReturnTrials==true

% ---------------- options ----------------
p = inputParser;
p.addParameter('ReturnTrials', false, @(x) islogical(x) || ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

% ---------------- basics ----------------
[P, K] = size(beta); %#ok<NASGU>
[N, nW] = size(toneOn_Go_tbyb);
assert(isvector(winCtrs) && numel(winCtrs)==nW, 'winCtrs must be 1 x nW.');
assert(size(Xz,1) == N*nW, 'Xz must have N*nW rows.');

goIdx   = trI.goI(:);
nogoIdx = trI.nogoI(:);
assert(numel(goIdx)==N && numel(nogoIdx)==N, 'trI.goI/nogoI must be [N x 1].');

HestAll = cell(1,K);

% ---------------- precompute F for each event type ----------------
[~,~,F_on_go]   = convolve_events_with_basis(toneOn_Go_tbyb,    basisToneOn,  'toneOnGo',   'assertCentered', true);
[~,~,F_off_go]  = convolve_events_with_basis(toneOff_Go_tbyb,   basisToneOff, 'toneOffGo',  'assertCentered', true);
[~,~,F_on_ng]   = convolve_events_with_basis(toneOn_NoGo_tbyb,  basisToneOn,  'toneOnNoGo', 'assertCentered', true);
[~,~,F_off_ng]  = convolve_events_with_basis(toneOff_NoGo_tbyb, basisToneOff, 'toneOffNoGo','assertCentered', true);

% Event-specific contributions (still in z-units)
H_on_go   = buildH_for_all_motifs(F_on_go,   'toneOnGo',   beta, X_names, muX, sdX);
H_off_go  = buildH_for_all_motifs(F_off_go,  'toneOffGo',  beta, X_names, muX, sdX);
H_on_ng   = buildH_for_all_motifs(F_on_ng,   'toneOnNoGo', beta, X_names, muX, sdX);
H_off_ng  = buildH_for_all_motifs(F_off_ng,  'toneOffNoGo',beta, X_names, muX, sdX);

% Full prediction from Xz * beta (z-units)
H_full = build_full_from_X(Xz, N, nW, beta);   % [N x nW x K]

% ---------------- summarize per motif ----------------
for k = 1:K
    Hk_on_go   = squeeze(H_on_go(:,:,k));    % [N x nW]
    Hk_off_go  = squeeze(H_off_go(:,:,k));
    Hk_on_ng   = squeeze(H_on_ng(:,:,k));
    Hk_off_ng  = squeeze(H_off_ng(:,:,k));
    Hk_full    = squeeze(H_full(:,:,k));

    % Event-specific means/SEMs
    [m_on_go,  s_on_go,  T_on_go ] = meanSemMasked(Hk_on_go,  goIdx,   opt.ReturnTrials);
    [m_off_go, s_off_go, T_off_go] = meanSemMasked(Hk_off_go, goIdx,   opt.ReturnTrials);
    [m_on_ng,  s_on_ng,  T_on_ng ] = meanSemMasked(Hk_on_ng,  nogoIdx, opt.ReturnTrials);
    [m_off_ng, s_off_ng, T_off_ng] = meanSemMasked(Hk_off_ng, nogoIdx, opt.ReturnTrials);

    % Full prediction means/SEMs
    [m_full_go,  s_full_go,  T_full_go ] = meanSemMasked(Hk_full, goIdx,   opt.ReturnTrials);
    [m_full_ng,  s_full_ng,  T_full_ng ] = meanSemMasked(Hk_full, nogoIdx, opt.ReturnTrials);

    Hk = struct();
    Hk.motifIdx = k;
    Hk.mean = struct( ...
        'toneOnGo',    m_on_go, ...
        'toneOffGo',   m_off_go, ...
        'toneOnNoGo',  m_on_ng, ...
        'toneOffNoGo', m_off_ng, ...
        'fullGo',      m_full_go, ...
        'fullNoGo',    m_full_ng);
    Hk.sem = struct( ...
        'toneOnGo',    s_on_go, ...
        'toneOffGo',   s_off_go, ...
        'toneOnNoGo',  s_on_ng, ...
        'toneOffNoGo', s_off_ng, ...
        'fullGo',      s_full_go, ...
        'fullNoGo',    s_full_ng);
    Hk.time = winCtrs;

    if opt.ReturnTrials
        Hk.trial = struct( ...
            'toneOnGo',    T_on_go, ...
            'toneOffGo',   T_off_go, ...
            'toneOnNoGo',  T_on_ng, ...
            'toneOffNoGo', T_off_ng, ...
            'fullGo',      T_full_go, ...
            'fullNoGo',    T_full_ng);
    end

    HestAll{k} = Hk;
end

end

%% =======================================================================
function H_all = buildH_for_all_motifs(F, baseName, beta, X_names, muX, sdX)
    [N, nW, nB] = size(F);
    [~, K]      = size(beta);

    % Find β columns for this event block, based on baseName_rc##
    pat = ['^' regexpquote(baseName) '_rc(\d+)$'];
    idx = find(~cellfun('isempty', regexp(X_names, pat, 'once')));

    if isempty(idx)
        warning('motif_event_estimates: No β columns found for %s; returning zeros.', baseName);
        H_all = zeros(N, nW, K);
        return;
    end

    % Sort by rc number (##) to match basis order
    rcnum = nan(size(idx));
    for ii = 1:numel(idx)
        tok = regexp(X_names{idx(ii)}, '_rc(\d+)$', 'tokens','once');
        if ~isempty(tok)
            rcnum(ii) = str2double(tok{1});
        end
    end
    [~, ord] = sort(rcnum, 'ascend');
    idx = idx(ord);

    nB = size(F,3);
    if numel(idx) ~= nB
        error('motif_event_estimates: %s has %d bases, but %d β columns found.', ...
              baseName, nB, numel(idx));
    end

    % Reshape F and standardize using stored μ/σ (same transform as Xz)
    F2 = reshape(F, [], nB);      % [N*nW x nB]
    muBlock = muX(idx);           % 1 x nB
    sdBlock = sdX(idx);           % 1 x nB
    sdBlock(~isfinite(sdBlock) | sdBlock < 1e-3) = 1e-3;

    F2z = (F2 - muBlock) ./ sdBlock;
    F2z(~isfinite(F2z)) = 0;

    H_all = zeros(N, nW, K);
    for k = 1:K
        w    = beta(idx, k);      % [nB x 1]
        Hvec = F2z * w;           % [N*nW x 1], z-units
        H_all(:,:,k) = reshape(Hvec, N, nW);
    end
end

function H_full = build_full_from_X(Xz, N, nW, beta)
% BUILD_FULL_FROM_X  Full GLM prediction Xz*beta reshaped to [N x nW x K].
    [M, ~] = size(Xz);
    [~, K] = size(beta);
    assert(M == N*nW, 'Xz rows must equal N*nW.');

    H_full = zeros(N, nW, K);
    for k = 1:K
        yhat = Xz * beta(:,k);          % [M x 1], z-units
        H_full(:,:,k) = reshape(yhat, N, nW);
    end
end

function [m, s, Trials] = meanSemMasked(Hk, mask, returnTrials)
    Hsel = Hk(mask, :);
    if isempty(Hsel)
        m = zeros(1, size(Hk,2));
        s = zeros(1, size(Hk,2));
        Trials = [];
    else
        m = mean(Hsel, 1, 'omitnan');
        s = std(Hsel, 0, 1, 'omitnan') ./ sqrt(size(Hsel,1));
        if returnTrials
            Trials = Hsel;
        else
            Trials = [];
        end
    end
end

function s = regexpquote(s)
    s = regexprep(s, '([.^$*+?{}\[\]\\|()])', '\\$1');
end
