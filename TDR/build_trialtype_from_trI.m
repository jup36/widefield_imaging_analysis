function TT = build_trialtype_from_trI(trI, varargin)
%BUILD_TRIALTYPE_FROM_TRI  Create categorical labels and design matrices.
%
% TT = BUILD_TRIALTYPE_FROM_TRI(trI, 'Coding', 'reference', 'Ref', 'CR')
%   trI fields expected (logical, N×1 or 1×N): hitI, missI, crI, faI
%
% OUTPUT (struct TT)
%   .labels        N×1 categorical with levels: Hit, Miss, CR, FA
%   .onehot        N×4 double, columns in order [Hit Miss CR FA] (0/1)
%   .X_ref         N×(4-1) double, reference coding (drops Ref)
%   .refLevel      char/string, the dropped reference (default 'CR')
%   .X_effect      N×(4-1) double, effect coding (sum-to-zero)
%   .names         1×4 string array of column names
%   .N             scalar, #trials
%
% OPTIONS
%   'Coding'   'reference'|'effect'|{'both'}  (only affects .X_* presence)
%   'Ref'      reference level for reference coding (default 'CR')
%
% Junchol-compatible helper, 2025.

p = inputParser;
p.addParameter('Coding','both', @(s) any(strcmpi(s,{'reference','effect','both'})));
p.addParameter('Ref','CR', @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

% ---- gather & sanitize
names = ["Hit","Miss","CR","FA"];
fns   = ["hitI","missI","crI","faI"];
assert(all(isfield(trI, cellstr(fns))), 'trI must contain hitI, missI, crI, faI.');

% coerce to column logicals
toCol = @(x) logical(x(:));
H = toCol(trI.hitI);
M = toCol(trI.missI);
C = toCol(trI.crI);
F = toCol(trI.faI);

N = numel(H);

% exclusivity / exhaustiveness checks
sumRows = double(H)+double(M)+double(C)+double(F);
assert(all(sumRows==1), 'Each trial must have exactly one of Hit/Miss/CR/FA true.');

% labels
lbl = strings(N,1);
lbl(H) = "Hit"; lbl(M) = "Miss"; lbl(C) = "CR"; lbl(F) = "FA";
labels = categorical(lbl, names);   % fixed order

% one-hot (N×4)
onehot = [double(H) double(M) double(C) double(F)]; % columns: Hit Miss CR FA

% reference coding (drop one column)
refLevel = string(opt.Ref);
assert(any(names==refLevel), 'Ref must be one of: %s', strjoin(names, ', '));
keepMask = names~=refLevel;
X_ref = onehot(:, keepMask);
refNames = names(keepMask);

% effect coding (sum-to-zero). Map levels to contrasts vs. grand mean.
% For K levels, produce K-1 columns; last level is -sum of others.
K = numel(names);
Cmat = eye(K) - (1/K);       % center columns (K×K)
Cmat(:,end) = [];            % keep first K-1 columns
% reorder rows to [Hit Miss CR FA]
X_effect = onehot * Cmat;    % N×(K-1)
effNames = names(1:end-1);

% assemble output
TT = struct();
TT.labels   = labels;
TT.onehot   = onehot;
TT.N        = N;
TT.names    = names;

switch lower(opt.Coding)
    case 'reference'
        TT.X_ref    = X_ref;
        TT.refLevel = refLevel;
    case 'effect'
        TT.X_effect = X_effect;
    otherwise  % both
        TT.X_ref    = X_ref;    TT.refLevel = refLevel;
        TT.X_effect = X_effect;
end

% (optional) convenience masks
TT.isGo    = toCol(isfield(trI,'goI')    && ~isempty(trI.goI)    .* trI.goI)    | H | M; % fallback
TT.isNoGo  = toCol(isfield(trI,'nogoI')  && ~isempty(trI.nogoI)  .* trI.nogoI)  | C | F;
TT.isRwd   = isfield(trI,'waterI')   && ~isempty(trI.waterI)   && logical(trI.waterI(:));
TT.isPun   = isfield(trI,'airpuffI') && ~isempty(trI.airpuffI) && logical(trI.airpuffI(:));
TT.isLicked= isfield(trI,'lickI')    && ~isempty(trI.lickI)    && logical(trI.lickI(:));

end
