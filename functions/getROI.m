%   This function outputs coordinates of the window from which to get the
%   aspiration depth
%
%
%   newframes are video frames read in previously
%   manualCorner = 0 if automatic corner detection is to be used
%                = 1 if ROI needs to be selected manually (could happen
%                    if there is a problem with edge detection)

%   coords = [xmin xmax ymin ymax]
%   ROIframes = newframes cropped to ROI size
%


function [coords ROIframes] = getROI(newframes, manualCorner, cannyThresh)

if nargin < 3
    cannyThresh = .45;
end

if manualCorner
    
    figure, imshow(newframes(:,:,1));
    title('draw region with top corner of pipette');
    coord = getrect;
    % coords = [xmin xmax ymin ymax]
    coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
    ROIframes = newframes(coords(3):coords(4), coords(1):coords(2), :);
    
else
    
%     figure, imshow(newframes(:,:,1))
    BW = edge(newframes(:,:,1), 'canny', cannyThresh);
    figure, imshow(BW);
    title('draw region with top corner of pipette');
    coord = getrect;
    % coords = [xmin xmax ymin ymax]
    coords =[coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)];
    
    % given rectangle with top corner of pipette, extract coordinates of edge
    % line and return "corner" point
    [xout yout] = returnCorner(BW, coords);
    close all;
    fig = figure;
    imshow(BW);
    truesize(fig, 2*size(newframes(:,:,1)));
    title('choose bottom right bound of ROI');
    [xb yb] = ginput(1);
    ROIframes = newframes(round(yout:yb), round(xout:xb), :);
    
end