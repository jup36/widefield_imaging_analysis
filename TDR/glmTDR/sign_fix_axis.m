function v = sign_fix_axis(v, mode)
mode = lower(string(mode));
switch mode
    case "maxabs"
        [~,idx] = max(abs(v));
        if v(idx) < 0, v = -v; end
    case "none"
        % no-op
    otherwise
        % no-op
end
end