function concatenateInitialFrames(folderPath)
    % This function takes a folder path as input and concatenates the first
    % 10 .png images named 'frame_1.png' to 'frame_10.png' horizontally.
    
    % Initialize an empty array to store the images
    concatenatedImage = [];
    
    % Loop through the first 10 frames
    for i = 1:10
        % Construct the filename
        fileName = fullfile(folderPath, sprintf('frame_%d.png', i));
        
        % Check if the file exists
        if exist(fileName, 'file')
            % Read the image
            img = imread(fileName);
            
            % Concatenate the image with the previous images
            concatenatedImage = [concatenatedImage, img];
        else
            warning('File %s does not exist. Skipping.', fileName);
        end
    end
    
    % If the concatenatedImage is not empty, save the image
    if ~isempty(concatenatedImage)
        % Construct the output filename
        outputFileName = fullfile(folderPath, 'initialFramesConcat.png');
        
        % Write the concatenated image to file
        imwrite(concatenatedImage, outputFileName);
        
        fprintf('Concatenated image saved as %s\n', outputFileName);
    else
        error('No images were found or concatenated.');
    end
end
