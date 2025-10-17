function TT = build_trial_masks(trI)
%BUILD_TRIAL_MASKS  Make trial-type labels and convenience masks from trI.
% Expects trI to contain logical vectors for some/all of:
%   hitI, missI, crI, faI, goI, nogoI, waterI, airpuffI, lickI
%
% Returns struct TT with:
%   .N             number of trials
%   .Hit/.Miss/.CR/.FA  (N×1 logical)
%   .isGo/.isNoGo/.isRwd/.isPun/.isLicked (N×1 logical)
%   .trialType     (N×1 categorical: 'Hit','Miss','CR','FA')
%   .idx           struct of linear indices for each class

% ----- resolve N (prefer hitI if present) -----
cand = {'hitI','missI','crI','faI','goI','nogoI','waterI','airpuffI','lickI'};
N = [];
for k = 1:numel(cand)
    if isfield(trI, cand{k}) && ~isempty(trI.(cand{k}))
        N = numel(trI.(cand{k}));
        break
    end
end
assert(~isempty(N), 'build_trial_masks:CannotInferN', ...
    'Cannot infer #trials from trI; no usable fields found.');

% ----- helper to fetch a column logical (or all-false if missing) -----
getMask = @(f) (isfield(trI,f) && ~isempty(trI.(f))) .* logical(trI.(f)(:));
toCol   = @(x) logical(x(:));

% ----- primary 4-class masks -----
H = getMask('hitI');     if isempty(H),     H = false(N,1); end
M = getMask('missI');    if isempty(M),     M = false(N,1); end
C = getMask('crI');      if isempty(C),     C = false(N,1); end
F = getMask('faI');      if isempty(F),     F = false(N,1); end

% ensure correct length
H = padOrTrim(H,N); M = padOrTrim(M,N); C = padOrTrim(C,N); F = padOrTrim(F,N);

% (optional) sanity: check exclusivity
overlap = (H&M) | (H&C) | (H&F) | (M&C) | (M&F) | (C&F);
if any(overlap)
    warning('build_trial_masks:Overlap', ...
        'Some trials are assigned to multiple classes; fixing by precedence Hit>Miss>CR>FA.');
    % enforce precedence Hit > Miss > CR > FA
    M(overlap & M & H) = false; C(H) = false; F(H) = false;
    C(M) = false; F(M) = false;
    F(C) = false;
end

% ----- convenience masks (safe fallbacks) -----
goI    = getMask('goI');     if isempty(goI),    goI    = false(N,1); end
nogoI  = getMask('nogoI');   if isempty(nogoI),  nogoI  = false(N,1); end
waterI = getMask('waterI');  if isempty(waterI), waterI = false(N,1); end
puffI  = getMask('airpuffI');if isempty(puffI),  puffI  = false(N,1); end
lickI  = getMask('lickI');   if isempty(lickI),  lickI  = false(N,1); end

% build composite masks
isGo     = goI   | H | M;
isNoGo   = nogoI | C | F;
isRwd    = waterI;
isPun    = puffI;
isLicked = lickI;

% ----- categorical labels -----
labels = repmat("NA", N, 1);
labels(F) = "FA";
labels(C) = "CR";
labels(M) = "Miss";
labels(H) = "Hit";
trialType = categorical(labels, {'Hit','Miss','CR','FA','NA'});

% ----- pack output -----
TT = struct();
TT.N = N;
TT.Hit  = H;
TT.Miss = M;
TT.CR   = C;
TT.FA   = F;

TT.isGo     = toCol(isGo);
TT.isNoGo   = toCol(isNoGo);
TT.isRwd    = toCol(isRwd);
TT.isPun    = toCol(isPun);
TT.isLicked = toCol(isLicked);

TT.trialType = trialType;
TT.idx = struct('Hit', find(H), 'Miss', find(M), 'CR', find(C), 'FA', find(F));
end

function x = padOrTrim(x, N)
% make x logical Nx1 (pad with false or trim)
x = logical(x(:));
if numel(x) < N
    x(end+1:N,1) = false;
elseif numel(x) > N
    x = x(1:N);
end
end
