%% Helper: make baseline save header
function fileheader_save = localMakeBaselineSaveHeader(fileheader)
% Convert task/baseline naming to a compact baseline-safe output label.
%
% Examples:
%   m1045_122424_task_day4-8_img      -> m1045_122424_base_day4-8_img
%   m1045_122424_baseline_day4-8_img  -> m1045_122424_base_day4-8_img
%   m1045_122424_baseline_img         -> m1045_122424_base_img

fileheader_save = char(fileheader);

fileheader_save = regexprep(fileheader_save, '_task_', '_base_');
fileheader_save = regexprep(fileheader_save, '_baseline_', '_base_');

% Fallbacks for cases like "..._task" or "..._baseline" without following underscore
fileheader_save = regexprep(fileheader_save, '_task$', '_base');
fileheader_save = regexprep(fileheader_save, '_baseline$', '_base');

end