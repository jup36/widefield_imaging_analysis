function RefitCuratedStaticBasisMotifs_ArrayTask_JP(manifest_fn)
% RefitCuratedStaticBasisMotifs_ArrayTask_JP
%
% SLURM array wrapper for static basis motif refitting.
% Each array task uses SLURM_ARRAY_TASK_ID as the chunk index.

fprintf('\n============================================================\n');
fprintf('RefitCuratedStaticBasisMotifs_ArrayTask_JP\n');
fprintf('Manifest: %s\n', manifest_fn);
fprintf('============================================================\n');

if exist(manifest_fn, 'file') ~= 2
    error('Array manifest not found: %s', manifest_fn);
end

S = load(manifest_fn);

requiredFields = { ...
    'file_processed_bucket', ...
    'basis_dir_bucket', ...
    'save_dir_bucket', ...
    'parameter_class', ...
    'nChunks'};

for i = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{i})
        error('Manifest missing required field: %s', requiredFields{i});
    end
end

taskIDstr = getenv('SLURM_ARRAY_TASK_ID');

if isempty(taskIDstr)
    error('SLURM_ARRAY_TASK_ID is empty. This must be run as a SLURM array task.');
end

chunk = str2double(taskIDstr);

if isnan(chunk) || chunk < 1 || chunk > S.nChunks
    error('Invalid SLURM_ARRAY_TASK_ID=%s for nChunks=%d.', taskIDstr, S.nChunks);
end

fprintf('\nArray task ID / chunk: %d/%d\n', chunk, S.nChunks);
fprintf('file_processed: %s\n', S.file_processed_bucket);
fprintf('basis_dir: %s\n', S.basis_dir_bucket);
fprintf('save_dir: %s\n', S.save_dir_bucket);
fprintf('parameter_class: %s\n', S.parameter_class);

RefitCuratedStaticBasisMotifs_JP( ...
    S.file_processed_bucket, ...
    S.basis_dir_bucket, ...
    chunk, ...
    S.parameter_class, ...
    S.save_dir_bucket);

fprintf('\nCompleted static refit array task for chunk %d/%d\n', chunk, S.nChunks);

end