function out = setNegNaN(in)
    in(in < 0) = NaN;
    out = in;
end