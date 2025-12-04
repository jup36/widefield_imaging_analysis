function kc = center_kernel_linear(k, idx0, centerIdx)
%CENTER_KERNEL_LINEAR  Shift kernel so sample at lag==0 lands at centerIdx.
% Zero-pad instead of wrapping to keep positive-lag energy on the right.
    L = numel(k);
    s = centerIdx - idx0;   % desired shift in samples (+ = right)
    if s > 0
        % shift right by s: prepend s zeros, drop last s samples
        kc = [zeros(s,1); k(1:L-s)];
    elseif s < 0
        % shift left by |s|: drop first |s| samples, append zeros
        s = -s;
        kc = [k(1+s:L); zeros(s,1)];
    else
        kc = k;
    end
end