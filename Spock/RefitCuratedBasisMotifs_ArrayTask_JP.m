function RefitCuratedBasisMotifs_ArrayTask_JP(manifest_fn)
% RefitCuratedBasisMotifs_ArrayTask_JP
%
% SLURM array wrapper for curated spatiotemporal basis motif refitting.
% Each array task uses SLURM_ARRAY_TASK_ID as the chunk index.
%
% This calls:
%   RefitCuratedBasisMotifs_JP

fprintf('\n============================================================\n');
fprintf('RefitCuratedBasisMotifs_ArrayTask_JP\n');
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
    'nChunks', ...
    'dateStr'};

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
fprintf('dateStr: %s\n', S.dateStr);

if isfield(S, 'H_init_method')
    H_init_method = S.H_init_method;
else
    H_init_method = 'projection';
end

if isfield(S, 'H_init_prctile_cap')
    H_init_prctile_cap = S.H_init_prctile_cap;
else
    H_init_prctile_cap = 99.5;
end

fprintf('H_init_method: %s\n', H_init_method);
fprintf('H_init_prctile_cap: %.2f\n', H_init_prctile_cap);

RefitCuratedBasisMotifs_JP( ...
    S.file_processed_bucket, ...
    S.basis_dir_bucket, ...
    chunk, ...
    S.parameter_class, ...
    S.save_dir_bucket, ...
    'dateStr', S.dateStr, ...
    'H_init_method', H_init_method, ...
    'H_init_prctile_cap', H_init_prctile_cap);

fprintf('\nCompleted curated refit array task for chunk %d/%d\n', chunk, S.nChunks);

end