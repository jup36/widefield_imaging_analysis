function generate_motif_gifs(data_filepath, output_dir, varargin)
% Generates and saves a GIF for each spatiotemporal motif found in W_basis.
%
% Args:
%   data_filepath (str): Path to the .mat file containing 
%                        the 'W_basis' array.
%   output_dir (str): Directory to save the generated GIF files.
%
% Optional Name-Value Args:
%   'cmap_name' (str): Name of the MATLAB colormap (e.g., 'viridis', 'parula', 'jet').
%                      Defaults to 'viridis'.
%   'vmin' (double): The minimum value for the color axis. Defaults to [].
%   'vmax' (double): The maximum value for the color axis. Defaults to [].
%   'frame_duration' (double): Duration of each frame in the GIF, in seconds.
%                              Defaults to 0.2 (200ms).

    % --- 1. Setup Input Parser for Optional Arguments ---
    p = inputParser;
    
    % Define default values
    default_cmap = 'viridis';
    default_vmin = [];
    default_vmax = [];
    default_duration = 0.2; % 200ms
    
    % Add optional arguments
    addParameter(p, 'cmap_name', default_cmap, @ischar);
    addParameter(p, 'vmin', default_vmin, @isnumeric);
    addParameter(p, 'vmax', default_vmax, @isnumeric);
    addParameter(p, 'frame_duration', default_duration, @isnumeric);
    
    % Parse the inputs
    parse(p, varargin{:});
    
    % Get results from the parser
    cmap_name = p.Results.cmap_name;
    vmin = p.Results.vmin;
    vmax = p.Results.vmax;
    frame_duration = p.Results.frame_duration;

    % --- 2. Setup and Load Data ---
    fprintf('Starting GIF generation...\n');
    fprintf('Loading data from: %s\n', data_filepath);
    
    % Load the .mat file
    try
        data_struct = load(data_filepath);
    catch ME
        fprintf('Error loading .mat file: %s\n', ME.message);
        return;
    end
    
    % Check for 'W_basis' variable
    if ~isfield(data_struct, 'W_basis')
        fprintf('Error: Could not find "W_basis" variable in %s.\n', data_filepath);
        fprintf('Found variables: %s\n', strjoin(fieldnames(data_struct), ', '));
        return;
    end
    W_basis = data_struct.W_basis;
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
       mkdir(output_dir);
       fprintf('Created output directory: %s\n', output_dir);
    end
    
    % Get the colormap matrix (e.g., 256x3 double)
    try
        cmap = colormap(cmap_name);
    catch ME
        fprintf('Warning: Colormap ''%s'' not found. Defaulting to ''parula''.\n', cmap_name);
        cmap = colormap('parula');
    end
    
    % --- 3. Validate and Reshape Data ---
    try
        [num_pixels, num_motifs, num_frames] = size(W_basis);
    catch ME
        fprintf('Error: W_basis does not have 3 dimensions. Shape is [%s]\n', ...
                num2str(size(W_basis)));
        return;
    end

    if num_pixels ~= 4096
        fprintf('Warning: Expected first dimension to be 4096 (64x64), but got %d.\n', ...
                num_pixels);
    end
    
    try
        % Reshape from (4096, 30, 10) -> (64, 64, 30, 10)
        reshaped_W = reshape(W_basis, 64, 64, num_motifs, num_frames);
        fprintf('Successfully loaded and reshaped data: (64, 64, %d motifs, %d frames)\n', ...
                num_motifs, num_frames);
    catch ME
        fprintf('Error: Could not reshape array. %s\n', ME.message);
        return;
    end

    % --- 4. Main Loop: Generate GIF per Motif ---
    for k = 1:num_motifs
        fprintf('  Processing motif %d/%d...\n', k, num_motifs);
        
        % Get all frames for this motif: shape (64, 64, 10)
        % FIX: Use squeeze() to remove the singleton dimension (size 1)
        % This changes motif_data from (64, 64, 1, 10) to (64, 64, 10)
        motif_data = squeeze(reshaped_W(:, :, k, :));
        
        % Determine color limits for this motif
        if isempty(vmin) || isempty(vmax)
            % Auto-calculate limits from this motif's data
            current_vmin = min(motif_data(:));
            current_vmax = max(motif_data(:));
        else
            % Use user-provided fixed limits
            current_vmin = vmin;
            current_vmax = vmax;
        end
        
        % Handle edge case where data is flat (min == max)
        if current_vmin == current_vmax
            current_vmin = current_vmin - 1e-9; % Add tiny offset
            current_vmax = current_vmax + 1e-9;
        end
        
        clim = [current_vmin, current_vmax];
        
        % Define the output filename for this GIF
        output_filename = fullfile(output_dir, sprintf('motif_%02d.gif', k));
        
        % --- 5. Apply Colormap and Save Each Frame ---
        for i = 1:num_frames
            % Get the (64, 64) frame
            % This will now correctly grab the i-th 2D frame
            frame = motif_data(:, :, i);
            
            % Convert the double-precision frame to a scaled uint8 indexed image
            % mat2gray scales data to [0, 1] based on clim
            % gray2ind converts [0, 1] to [1, N] indices for the colormap
            % N = size(cmap, 1), which is typically 256
            
            % Manually scale and clip to be robust
            frame_scaled = (frame - clim(1)) / (clim(2) - clim(1));
            frame_scaled(frame_scaled < 0) = 0; % Clip min
            frame_scaled(frame_scaled > 1) = 1; % Clip max
            
            % Convert to colormap index.
            % We want indices from 1 to size(cmap, 1)
            N_colors = size(cmap, 1);
            frame_indexed = uint8(floor(frame_scaled * (N_colors - 1)) + 1);

            % Write to GIF
            if i == 1
                % First frame: create the file
                imwrite(frame_indexed, cmap, output_filename, 'gif', ...
                        'Loopcount', inf, 'DelayTime', frame_duration);
            else
                % Subsequent frames: append to the file
                imwrite(frame_indexed, cmap, output_filename, 'gif', ...
                        'WriteMode', 'append', 'DelayTime', frame_duration);
            end
        end % end frame loop
        
        fprintf('    Saved %s\n', output_filename);
        
    end % end motif loop
    
    fprintf('GIF generation complete.\n');
end

% --- EXAMPLE USAGE ---
%
% % 1. Create dummy data and save it
% disp('Creating and saving dummy W_basis.mat file...');
% W_basis = rand(4096, 30, 10); % 4096 pixels, 30 motifs, 10 frames
% save('W_basis.mat', 'W_basis');
% 
% % 2. Define output directory
% output_dir = 'motif_gifs';
% 
% % 3. Call the function with default options
% generate_motif_gifs('W_basis.mat', output_dir);
% 
% % 4. Call the function with custom options
% % Using 'jet' colormap, fixed color axis, and faster frame rate
% % generate_motif_gifs('W_basis.mat', output_dir, ...
% %     'cmap_name', 'jet', ...
% %     'vmin', 0.2, ...
% %     'vmax', 0.8, ...
% %     'frame_duration', 0.1);
%