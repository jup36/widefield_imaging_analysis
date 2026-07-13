function y = local_axisRMS_overTime(A)
%LOCAL_AXISRMS_OVERTIME
% A is [nAxes x nTime].
% Returns one [1 x nTime] trajectory:
%
%   y(t) = sqrt(mean(A(:,t).^2, 'omitnan'))
%
% This treats selected axes as a subspace/vector and summarizes the
% magnitude of the slope vector at each time bin.

    if isempty(A)
        y = [];
        return;
    end

    A = double(A);

    if isvector(A)
        % If only one axis survives, preserve as 1 x T.
        y = abs(A(:)');
    else
        y = sqrt(mean(A.^2, 1, 'omitnan'));
    end
end