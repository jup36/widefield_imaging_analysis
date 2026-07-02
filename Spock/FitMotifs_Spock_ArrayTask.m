function FitMotifs_Spock_ArrayTask(manifest_fn)
% FitMotifs_Spock_ArrayTask
%
% Wrapper for SLURM array motif fitting.
% Each array task reads SLURM_ARRAY_TASK_ID and runs one chunk.

fprintf('\n============================================================\n');
fprintf('FitMotifs_Spock_ArrayTask\n');
fprintf('Manifest: %s\n', manifest_fn);
fprintf('============================================================\n');

if exist(manifest_fn, 'file') ~= 2
    error('Array manifest not found: %s', manifest_fn);
end

S = load(manifest_fn);

requiredFields = {'file_processed_bucket', 'save_fn_bucket', 'parameter_class', 'nChunks'};
for i = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{i})
        error('Manifest missing required field: %s', requiredFields{i});
    end
end

taskIDstr = getenv('SLURM_ARRAY_TASK_ID');

if isempty(taskIDstr)
    error('SLURM_ARRAY_TASK_ID is empty. This function must be run as a SLURM array task.');
end

chunk = str2double(taskIDstr);

if isnan(chunk) || chunk < 1 || chunk > S.nChunks
    error('Invalid SLURM_ARRAY_TASK_ID=%s for nChunks=%d.', taskIDstr, S.nChunks);
end

file_processed = S.file_processed_bucket;
save_fn = S.save_fn_bucket{chunk};
parameter_class = S.parameter_class;

fprintf('\nArray task ID: %d/%d\n', chunk, S.nChunks);
fprintf('file_processed: %s\n', file_processed);
fprintf('save_fn: %s\n', save_fn);
fprintf('parameter_class: %s\n', parameter_class);

FitMotifs_Spock(file_processed, save_fn, chunk, parameter_class);

fprintf('\nCompleted FitMotifs_Spock_ArrayTask for chunk %d/%d\n', chunk, S.nChunks);

end