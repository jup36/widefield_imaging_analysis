function pethStr = convert_pethTime_string(pethTime)
    % Convert each number to a string format with 'p' for decimal and 'n' for negative values
    pethStrArray = arrayfun(@(x) format_number(x), pethTime, 'UniformOutput', false);
    
    % Concatenate with '_'
    pethStr = strjoin(pethStrArray, '_');
end

function formattedStr = format_number(num)
    % Convert number to string with one decimal place
    if num < 0
        formattedStr = sprintf('n%d', abs(floor(num))); % Negative prefix
    else
        formattedStr = sprintf('%d', floor(num));
    end
    
    % Append decimal part
    decimalPart = abs(num - floor(num));
    if decimalPart > 0
        formattedStr = sprintf('%sp%d', formattedStr, round(decimalPart * 10));
    else
        formattedStr = sprintf('%sp0', formattedStr);
    end
end