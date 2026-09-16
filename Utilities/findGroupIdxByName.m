function idx = findGroupIdxByName(groupsCell, targetName)
% Find the row index into cvR2_unique_grp / cvR2_upper_grp corresponding
% to a named group, by matching opts.Groups{g}.name. Returns [] if not
% found (e.g. that group had zero columns and was filtered out for that
% session -- shouldn't happen for GoToneOn/NoGoToneOn since they're always
% built, but checked defensively rather than assumed).
idx = [];
for g = 1:numel(groupsCell)
    if isfield(groupsCell{g}, 'name') && strcmp(groupsCell{g}.name, targetName)
        idx = g;
        return;
    end
end
end