% Inputs: path to where video files are stored, and current file name
%
% This function takes a time-lapse video of an embryo and outputs a vector
% of image processing params for input to an SVM
%
% Livia Zarnescu
% 3-25-15
%

function paramOutVector = extractEmbryoParams(videoPath, currFileName)


% read in video, keep only first color channel
obj = VideoReader([videoPath currFileName]);
frameRate = obj.FrameRate;
frames = read(obj);
grayframes = squeeze(frames(:,:,1,:));
s = size(grayframes);

% start extracting params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. detect approximately where in the cell the embryo starts out, and
% track its trajectory over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. define an ROI around embryo trajectory, and extract some params from
% that ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










