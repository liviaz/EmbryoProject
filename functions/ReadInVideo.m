% read in a video file and convert it to grayscale

function [framesOut, frameRate] = ReadInVideo(movpath, vidSecs, startFrame, flip)

% read in video
obj = VideoReader(movpath);
frameRate = obj.FrameRate;
totalFrames = round(frameRate*vidSecs);
obj.CurrentTime = startFrame / obj.FrameRate;
lastFrame = min(obj.FrameRate*obj.Duration - 1, startFrame + totalFrames);
numFrames = lastFrame - startFrame + 1;


% convert to grayscale and double format
framesOut = zeros(obj.Height, obj.Width, numFrames);

for i = startFrame:lastFrame
    
    currFrame = readFrame(obj);
    s = size(currFrame);
    
    if length(s) > 2
        % take each frame and make it grayscale
        if flip
            framesOut(:,:,i-startFrame+1) = rot90(double(rgb2gray(currFrame)),2)/255;
        else
            framesOut(:,:,i-startFrame+1) = double(rgb2gray(currFrame))/255;
        end
    else
        if flip
            framesOut(:,:,i-startFrame+1) = rot90(double(currFrame),2)/255;
        else
            framesOut(:,:,i-startFrame+1) = double(currFrame)/255;
        end
    end
end





