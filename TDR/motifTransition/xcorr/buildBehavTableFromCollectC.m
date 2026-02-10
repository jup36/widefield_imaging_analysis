function behTbl = buildBehavTableFromCollectC(sessInfo, hitTcollectC, crTcollectC, dPrmTcollectC, varargin)
% buildBehavTableFromCollectC
%   Build a per-session behavioral table aligned to sessInfo.header.
%
%   behTbl = buildBehavTableFromCollectC(sessInfo, hitTcollectC, crTcollectC, dPrmTcollectC, ...)
%
% INPUTS
%   sessInfo      : table with variable 'header' (e.g., "m1045_122424") in row order of Y
%   hitTcollectC  : Nmouse x 2 cell: col1=headers (cellstr), col2=hitRate (double)
%   crTcollectC   : Nmouse x 2 cell: col1=headers (cellstr), col2=crRate (double)
%   dPrmTcollectC : Nmouse x 2 cell: col1=headers (cellstr), col2=dPrime (double) [optional; can be []]
%
% NAME–VALUE
%   'clipEps'     : clip rates to [eps, 1-eps] before norminv (default 1e-3)
%   'computeDprimeFromRates' : if true, compute d' from hit/FA even if dPrmTcollectC provided (default false)
%
% OUTPUT
%   behTbl : table with variables
%       header, hitRate, crRate, faRate, dPrime, biasC
%
% NOTES
%   biasC is SDT criterion: c = -0.5*(Z(H) + Z(FA))
%

p = inputParser;
p.addParameter('clipEps', 1e-3, @(x) isnumeric(x) && isscalar(x) && x>0 && x<0.1);
p.addParameter('computeDprimeFromRates', false, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(istable(sessInfo) && ismember('header', sessInfo.Properties.VariableNames), ...
    'sessInfo must be a table containing variable ''header''.');

hdrSess = string(sessInfo.header(:));

% ---------- flatten collectC into maps ----------
hitMap = local_collectC_to_map(hitTcollectC);  % containers.Map header->value
crMap  = local_collectC_to_map(crTcollectC);

dMap = [];
if ~isempty(dPrmTcollectC)
    dMap = local_collectC_to_map(dPrmTcollectC);
end

S = numel(hdrSess);
hit = nan(S,1);
cr  = nan(S,1);
dp  = nan(S,1);

for i = 1:S
    key = char(hdrSess(i));
    if isKey(hitMap, key), hit(i) = hitMap(key); end
    if isKey(crMap,  key), cr(i)  = crMap(key);  end
    if ~isempty(dMap) && isKey(dMap, key), dp(i) = dMap(key); end
end

fa = 1 - cr;

% ---------- compute bias + (optionally) d' from rates ----------
eps0 = opt.clipEps;

hitC = hit; faC = fa;
validRate = isfinite(hitC) & isfinite(faC);
hitC(validRate) = min(max(hitC(validRate), eps0), 1-eps0);
faC(validRate)  = min(max(faC(validRate),  eps0), 1-eps0);

zH  = nan(S,1); zFA = nan(S,1);
zH(validRate)  = norminv(hitC(validRate));
zFA(validRate) = norminv(faC(validRate));

biasC = nan(S,1);
biasC(validRate) = -0.5 * (zH(validRate) + zFA(validRate));

if opt.computeDprimeFromRates || isempty(dMap)
    dp2 = nan(S,1);
    dp2(validRate) = zH(validRate) - zFA(validRate);
    dp = dp2;
end

behTbl = table();
behTbl.header   = hdrSess;
behTbl.hitRate  = hit;
behTbl.crRate   = cr;
behTbl.faRate   = fa;
behTbl.dPrime   = dp;
behTbl.biasC    = biasC;

end

% ===================== local helper =====================
function mp = local_collectC_to_map(collectC)
% collectC: N x 2 cell, {headers, values}
mp = containers.Map('KeyType','char','ValueType','double');

if isempty(collectC), return; end

for r = 1:size(collectC,1)
    hdrs = collectC{r,1};
    vals = collectC{r,2};

    if isempty(hdrs) || isempty(vals), continue; end

    hdrs = string(hdrs(:));
    vals = double(vals(:));

    n = min(numel(hdrs), numel(vals));
    for i = 1:n
        key = char(strtrim(hdrs(i)));
        if strlength(key)==0, continue; end
        mp(key) = vals(i);
    end
end
end
