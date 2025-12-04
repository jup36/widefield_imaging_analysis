function [B, lagsSec, info] = make_rcos_basis_ortho(nBases, lagRangeSec, binSec, varargin)
%MAKE_RCOS_BASIS_ORTHO  Raised-cosine temporal basis.
% Supports 'linear' and 'log' constructions. For 'log', matches the
% log_cos recipe: centers evenly in log-time with a 2-gap late margin and
% support ±2*gap in log space. Optional centerZero pads so lag==0 is at
% conv anchor. Optional orthonorm after construction.
%
% [B, lagsSec, info] = make_rcos_basis_ortho(nBases, [lagMin lagMax], binSec, ...)

% ---------- parse ----------
p = inputParser;
p.addParameter('nonlin','linear',@(s) any(strcmpi(s,{'linear','log'})));
p.addParameter('c',0.45,@(x) isnumeric(x) && isscalar(x) && x>0);   % width scale (linear only)
p.addParameter('normCols',true,@(v) islogical(v) || ismember(v,[0 1]));
p.addParameter('orthonorm',false,@(v) islogical(v) || ismember(v,[0 1]));
p.addParameter('centersSec',[],@(v) isempty(v) || (isnumeric(v) && isvector(v)));
p.addParameter('title','',@(s) isstring(s) || ischar(s));
p.addParameter('centerZero',false,@(v) islogical(v) || ismember(v,[0 1]));
% log options (matching your log_cos)
p.addParameter('logBiasFrac',0.10,@(x) isnumeric(x) && isscalar(x) && x>0); % fraction of window
p.addParameter('uLogWidthScale',1,@(x) isnumeric(x) && isscalar(x) && x>0);   % width scale (log only)
%uLogWidthScale = 1 → baseline overlap (~50%, like your current plot)
%uLogWidthScale > 1 → more overlap (broader kernels)
%uLogWidthScale < 1 → less overlap (narrower, more separated kernels)
p.addParameter('normalizeLogCols',true,@(v) islogical(v) || ismember(v,[0 1]));
p.parse(varargin{:});
opt = p.Results;

lagMin = lagRangeSec(1);
lagMax = lagRangeSec(2);
if lagMax <= lagMin, error('lagMax must be greater than lagMin.'); end

% ---------- lag grid ----------
lagsSec = (lagMin:binSec:lagMax).';
L = numel(lagsSec);

% ensure 0 is exactly on the grid
tol = max(1e-12, binSec*1e-6);
idx0 = find(abs(lagsSec) <= tol, 1, 'first');
if isempty(idx0)
    [~, idx0] = min(abs(lagsSec));
    lagsSec(idx0) = 0;
end

% ---------- build B ----------
B = [];
ctrSec = [];

switch lower(opt.nonlin)

case 'linear'
    % centers
    if ~isempty(opt.centersSec)
        ctrSec = opt.centersSec(:).';
        if numel(ctrSec) ~= nBases
            error('centersSec must have nBases elements.');
        end
    else
        ctrSec = linspace(lagMin, lagMax, nBases);
    end

    % width (linear)
    if nBases == 1
        w = (lagMax - lagMin) / 2;
    else
        w = opt.c * (ctrSec(2) - ctrSec(1));   % ~50% overlap when c~0.45–0.6
    end
    if w <= 0, error('Computed width w <= 0; check inputs.'); end

    % bumps in linear time
    B = zeros(L, nBases);
    halfw = w/2;
    for k = 1:nBases
        x = (lagsSec - ctrSec(k)) / halfw;
        x = max(-pi, min(pi, x));
        bump = 0.5*(cos(x)+1);
        bump(abs((lagsSec - ctrSec(k))/halfw) > pi) = 0;
        B(:,k) = bump;
    end

case 'log'
    % EXACT log_cos behavior (causal window)
    if lagMin < 0
        error('nonlin="log" requires lagRangeSec >= 0. Use a separate basis for negative lags.');
    end

    % bias in window units (like your bias = bias_min% * diff(range))
    win = max(lagMax - lagMin, binSec);
    bias = max(opt.logBiasFrac * win, 10*eps);  % e.g., 0.1 → 10% of window

    % warp to u-space
    u    = log(lagsSec + bias);
    uMin = log(lagMin + bias);
    uMax = log(lagMax + bias);

    % gaps and peaks in u with a 2-gap late margin
    gap_u  = (uMax - uMin) / (nBases + 1);
    u_peaks = uMin : gap_u : (uMax - 2*gap_u);     % 1 x nBases
    ctrSec  = exp(u_peaks) - bias;                  % for info only

    % build bumps with support ±2*gap_u in u; cosine argument scaled by pi/(2*gap_u)
    B = zeros(L, nBases);
    for k = 1:nBases 
        a = (u - u_peaks(k)) * (pi / (2*gap_u * opt.uLogWidthScale));  % (log_func(x)-peak)*pi/(2*gap*scaling factor)
        a = max(-pi, min(pi, a));
        bump = 0.5*(cos(a)+1);
        % zero outside ±2*gap in u
        outside = abs(u - u_peaks(k)) > 2*gap_u * opt.uLogWidthScale;
        bump(outside) = 0;
        B(:,k) = bump;
    end

    % optional column normalization (matches your log_cos 'normalize')
    if opt.normalizeLogCols
        s = sum(B,1); s(s==0)=1; B = B ./ s;
    end

otherwise
    error('Unknown nonlin=%s', opt.nonlin);
end

% ---------- centerZero (pad, no wrap/no crop) ----------
leftPad = 0; rightPad = 0; centerIdx = ceil((L+1)/2);
if opt.centerZero
    % want centerIdx == idx0 + leftPad, final length odd
    leftPad  = max(0, L + 1 - 2*idx0);
    rightPad = (2*idx0 - 1 - L) + leftPad;
    if leftPad>0 || rightPad>0
        B = [zeros(leftPad, nBases); B; zeros(rightPad, nBases)];
        lagsSec = [ (lagsSec(1) - (leftPad:-1:1)'*binSec ); lagsSec; (lagsSec(end) + (1:rightPad)'*binSec ) ];
        % re-pin zero
        [~, idx0] = min(abs(lagsSec)); lagsSec(idx0)=0;
        L = numel(lagsSec);
        centerIdx = ceil((L+1)/2);
    end
end

% ---------- orthonorm / norm ----------
if opt.orthonorm && nBases > 1
    [Q,~] = qr(B,0);
    for k = 1:size(Q,2)   % deterministic sign
        col = Q(:,k); j = find(abs(col)>10*eps,1,'first');
        if ~isempty(j) && col(j)<0, Q(:,k) = -Q(:,k); end
    end
    B = Q;
elseif opt.normCols
    s = sqrt(sum(B.^2,1)); s(s==0)=1; B = B ./ s;
end

% ---------- optional plot ----------
if ~isempty(opt.title)
    figure('Color','w');
    plot(lagsSec, B, 'LineWidth',1.2);
    xlabel('Lag (s)'); ylabel('Basis amplitude'); grid on
    title(opt.title, 'Interpreter','none');
end

% ---------- info ----------
info = struct('centersSec',ctrSec, 'lagsSec',lagsSec, ...
              'nonlin',opt.nonlin, 'c',opt.c, ...
              'normCols',logical(opt.normCols), 'orthonorm',logical(opt.orthonorm), ...
              'centerZero',logical(opt.centerZero), ...
              'leftPad',leftPad,'rightPad',rightPad,'centerIdx',centerIdx,'idx0',idx0);
end
