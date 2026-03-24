function idx = axisNameContains(axisNames, keyword)
% axisNameContains
% Returns logical index of axis names containing the keyword.
%
% INPUT
%   axisNames : cell array of strings (e.g. rezProjStats.perMouse{...}.axes.names)
%   keyword   : string or char
%
% OUTPUT
%   idx       : logical index vector

    axisNames = string(axisNames);
    keyword   = string(keyword);

    idx = contains(axisNames, keyword);
end