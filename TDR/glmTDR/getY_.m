function [Y, N, K] = getY_(glmRez, whichY, nW)
whichY = lower(string(whichY));
switch whichY
    case "yz"
        assert(isfield(glmRez,'Yz') && ~isempty(glmRez.Yz), 'glmRez.Yz missing/empty.');
        Y = glmRez.Yz;
    case "ybig"
        assert(isfield(glmRez,'Ybig') && ~isempty(glmRez.Ybig), 'glmRez.Ybig missing/empty.');
        Y = glmRez.Ybig;
    otherwise
        error('Unknown ProjectWhichY: %s', whichY);
end

assert(isnumeric(Y) && ndims(Y)==2, 'Y must be [M x K] numeric.');
[M,K] = size(Y);
assert(rem(M,nW)==0, 'M=%d not divisible by nW=%d -> stacking assumption broken.', M, nW);
N = M/nW;
end
