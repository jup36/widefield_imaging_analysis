function [earlyJ, lateJ] = splitHalves(validJ, dropMiddleWhenOdd)
% Splits an ordered vector of valid session indices into first/second half.
% Odd count: middle element dropped (or given to the late half).
n = numel(validJ);
if n == 0
    earlyJ = []; lateJ = [];
    return;
end

if mod(n, 2) == 0
    h = n / 2;
    earlyJ = validJ(1:h);
    lateJ  = validJ(h+1:end);
else
    h = floor(n / 2);
    earlyJ = validJ(1:h);
    if dropMiddleWhenOdd
        lateJ = validJ(h+2:end);      % skip the middle session
    else
        lateJ = validJ(h+1:end);      % middle session joins the late half
    end
end
end