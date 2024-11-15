function h = plotFrameSequence(dff, frameIndices, varargin)

if nargin < 3
    climit = [-2 2]; 
elseif nargin == 3
    climit = varargin{1}; 
end

% Create a figure
h = figure;

% Get the dimensions of the input matrix
[nX, nY, ~] = size(dff);

% Extract the frames specified by frameIndices
selectedFrames = dff(:,:,frameIndices);

% Initialize an empty matrix to store concatenated frames
concatenatedFrames = zeros(nX, nY * length(frameIndices));

% Loop through the frames and concatenate them horizontally
for i = 1:length(frameIndices)
    concatenatedFrames(:, (i-1)*nY + 1 : i*nY) = selectedFrames(:,:,i);
end

% Display the concatenated frames
imagesc(concatenatedFrames);
%clim(climit); 
colormap('parula');
colorbar; % Optional: adds a color bar to indicate the value scale
pbaspect([length(frameIndices) 1 1]);

end
