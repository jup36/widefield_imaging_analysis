function T = localRowsToTable(rows, varNames)
% Build a table from cell rows. Handles empty row sets safely.

if isempty(rows)
    rows = cell(0, numel(varNames));
end

T = cell2table(rows, 'VariableNames', varNames);

end