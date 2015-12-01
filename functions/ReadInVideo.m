% read in a video file and convert it to grayscale

function [framesOut, frameRate] = ReadInVideo(movpath, vidSecs, startFrame, flip)

% read in video
obj = VideoReader(movpath);
frameRate = obj.FrameRate;
totalFrames = round(frameRate*vidSecs);
lastFrame = min(obj.NumberOfFrames - 1, startFrame + totalFrames);
frames = read(obj, [startFrame lastFrame]);

% convert frames to grayscale and double format
s = size(frames);
numFrames = s(4);
framesOut = zeros(size(frames,1), size(frames,2), numFrames);
for i = 1:numFrames
    if s(3) > 1
        % take each frame and make it grayscale
        if flip
            framesOut(:,:,i) = fliplr(double(rgb2gray(frames(:,:,:,i))))/255;
        else
            framesOut(:,:,i) = double(rgb2gray(frames(:,:,:,i)))/255;
        end
    else
        if flip
            framesOut(:,:,i) = fliplr(double(frames(:,:,i)))/255;
        else
            framesOut(:,:,i) = double(frames(:,:,i))/255;
        end
    end
end





