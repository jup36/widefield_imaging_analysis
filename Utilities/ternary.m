% ===== Local functions (must be after all script code) =====
function y = ternary(cond, a, b)
%TERNARY Simple inline-style conditional.
%   y = TERNARY(cond, a, b) returns a if cond is true, else b.
%   Works for strings, numerics, structs, cells, etc.
    if cond
        y = a;
    else
        y = b;
    end
end
