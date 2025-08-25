function v = getOpt(obj, name, defaultVal)
try
    if isstruct(obj) && isfield(obj, name)
        v = obj.(name);
    elseif isobject(obj) && isprop(obj, name)
        v = obj.(name);
    else
        v = defaultVal;
    end
catch
    v = defaultVal;
end
end