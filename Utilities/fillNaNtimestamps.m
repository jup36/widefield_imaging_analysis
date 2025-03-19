function filledTimestamps = fillNaNtimestamps(timestamps)
% FILLNANTIMESTAMPS - Fills NaN values in a timestamp array by identifying consecutive NaN chunks 
% and retrospectively filling them based on a strict 100 ms interval.
%
% Inputs:
%   timestamps - A numeric array with NaN values representing missing timestamps.
%
% Outputs:
%   filledTimestamps - The array with NaN values replaced based on a 100 ms interval.

    % Define the expected interval (100 ms = 0.1 seconds)
    expectedInterval = 0.1;  

    % Copy the original timestamps for modification
    filledTimestamps = timestamps;

    % Iterate through the timestamps
    idx = 1;
    while idx <= length(filledTimestamps)
        % If we encounter a NaN
        if isnan(filledTimestamps(idx))
            startNaN = idx; % Record the start of the NaN chunk
            
            % Find the first non-NaN value after the NaN chunk
            while idx <= length(filledTimestamps) && isnan(filledTimestamps(idx))
                idx = idx + 1;
            end
            
            if idx > length(filledTimestamps)
                % If we reached the end of the array without finding a valid timestamp, stop
                break;
            end
            
            % The first valid timestamp after the NaN chunk
            nextValidTime = filledTimestamps(idx);
            
            % Compute the missing timestamps by stepping backwards from the valid timestamp
            numMissing = idx - startNaN;
            filledTimestamps(startNaN:idx-1) = nextValidTime - (numMissing:-1:1) * expectedInterval;
        else
            % Move to the next index
            idx = idx + 1;
        end
    end
end
