function figHandle = imageFrameWithNaNsEdgeMap(frame, climits, edgeMap, edgeColor, edgeAlpha)
% climits = [-3 3];

% Step 1: Determine the rows and columns that contain non-zero values
[y, x] = find(edgeMap);

% Step 2: Find the bounding box
minRow = max(min(y) - 2, 1); % Ensure not going out of bounds
maxRow = min(max(y) + 2, size(frame, 1)); 
minCol = max(min(x) - 2, 1);
maxCol = min(max(x) + 2, size(frame, 2));

% Step 3: Crop the images
croppedFrame = frame(minRow:maxRow, minCol:maxCol);
croppedEdgeMap = edgeMap(minRow:maxRow, minCol:maxCol);

% Create a new figure and store its handle
figHandle = figure;

% Create an adjusted colormap with white at the beginning
cmap = [0, 0, 0; hot(256)];

% Display the cropped frame
imagesc(croppedFrame);

% Use the custom colormap
colormap(cmap);

% Set NaN values to a value outside the data range for visualization
clim(climits);

% Hold on to overlay the cropped edge map
hold on;

% Extract x and y coordinates of the edges from the croppedEdgeMap
[yCropped, xCropped] = find(croppedEdgeMap);

% Plot each edge using a tiny square patch to represent a point
for i = 1:length(xCropped)
    patch([xCropped(i)-0.5, xCropped(i)+0.5, xCropped(i)+0.5, xCropped(i)-0.5], [yCropped(i)-0.5, yCropped(i)-0.5, yCropped(i)+0.5, yCropped(i)+0.5], edgeColor, 'FaceAlpha', edgeAlpha, 'EdgeColor', 'none');
end

% Show the colorbar (optional)
colorbar;

hold off;

end
