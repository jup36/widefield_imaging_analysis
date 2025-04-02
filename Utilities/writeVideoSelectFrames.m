function writeVideoSelectFrames(filePathOrg, frameVec, filePathOut)
% writeVideoSelectFrames extracts selected frames from a video and writes them to a new video file.
%
% Inputs:
%   - filePathOrg: full path to the original video file
%   - frameVec: vector of frame indices to extract (e.g., [1 5 10 30])
%   - filePathOut: full path to save the new video
%
% Tip: Use .avi with 'Motion JPEG AVI' for safe intra-frame video writing.

    % Input validation
    if nargin ~= 3
        error('Usage: writeVideoSelectFrames(filePathOrg, frameVec, filePathOut)');
    end
    if ~isfile(filePathOrg)
        error('Original video not found: %s', filePathOrg);
    end
    if ~isvector(frameVec) || any(frameVec < 1) || any(mod(frameVec,1) ~= 0)
        error('frameVec must be a vector of positive integers.');
    end

    % Load video
    vidReader = VideoReader(filePathOrg);
    totalFrames = floor(vidReader.Duration * vidReader.FrameRate);

    frameVec = sort(unique(frameVec));
    if any(frameVec > totalFrames)
        error('Frame index exceeds total number of frames (%d)', totalFrames);
    end

    % Use safe codec: Motion JPEG AVI (each frame encoded independently)
    if endsWith(filePathOut, '.mp4')
        warning('Consider using .avi instead of .mp4 to avoid compression artifacts.');
    end
    vidWriter = VideoWriter(filePathOut, 'Motion JPEG AVI');
    vidWriter.FrameRate = vidReader.FrameRate;
    open(vidWriter);

    fprintf('Writing selected frames to: %s\n', filePathOut);

    % Initialize
    nextIdx = 1;
    nextFrameToWrite = frameVec(nextIdx);
    frameCounter = 0;

    % Loop through video sequentially
    while hasFrame(vidReader)
        frame = readFrame(vidReader);
        frameCounter = frameCounter + 1;

        if frameCounter == nextFrameToWrite
            writeVideo(vidWriter, frame);
            fprintf('  → Wrote frame %d (%d of %d)\n', ...
                frameCounter, nextIdx, numel(frameVec));

            nextIdx = nextIdx + 1;
            if nextIdx > numel(frameVec)
                break;  % Done writing all requested frames
            end
            nextFrameToWrite = frameVec(nextIdx);
        end
    end

    close(vidWriter);
    fprintf('✅ Done. Wrote %d frames to %s\n', numel(frameVec), filePathOut);
end
