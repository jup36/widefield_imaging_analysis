function trialMask = local_trialMask_fromGroup_(trId, groupName, N, policy)
policy = string(policy);
groupName = string(groupName);

trialMask = true(N,1);
want = "all";

if strcmpi(policy, "GoNoGoByGroup")
    if contains(lower(groupName), "nogo")
        want = "nogo";
    else
        want = "go";
    end
elseif strcmpi(policy, "CorrectOnlyByGroup")
    if contains(lower(groupName), "nogo")
        want = "cr";
    else
        want = "hit";
    end
end

if want == "all" || isempty(trId) || ~isstruct(trId)
    return;
end

switch want
    case "go",   trialMask = local_makeMask_(trId, 'goI',   N);
    case "nogo", trialMask = local_makeMask_(trId, 'nogoI', N);
    case "hit",  trialMask = local_makeMask_(trId, 'hitI',  N);
    case "cr",   trialMask = local_makeMask_(trId, 'crI',   N);
end

if ~any(trialMask)
    warning('local_trialMask_fromGroup_:EmptyMask', ...
        'Empty trial mask for %s (%s); using all trials', groupName, want);
    trialMask = true(N,1);
end
end