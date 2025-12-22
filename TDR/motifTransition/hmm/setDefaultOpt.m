function opt = setDefaultOpt(opt, field, val)
if ~isfield(opt, field) || isempty(opt.(field))
    opt.(field) = val;
end
end