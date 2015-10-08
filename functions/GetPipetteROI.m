% get ROI inside pipette
% do it based on template matching, not just corner detection anymore

function [ROIframesOut] = GetPipetteROI(frames, ...
    cannyThresh, extraFig, filePathRaw)

if isequal(filePathRaw(end), '\') || isequal(filePathRaw(end), '/')
    filePathRaw = filePathRaw(1:end-1);
end

if ~exist([filePathRaw '\pipRef.mat'], 'file')
    % save template of original 
    MakeRoiTemplate(filePathRaw, cannyThresh, extraFig);
end


% template matching from mat file
load([filePathRaw '\pipRef.mat']);

[~,I_NCC,~] = template_matching(ROIframe, frames(:,:,1));
[yBest,xBest] = find(I_NCC == max(I_NCC(:)));

yCrop = round(yBest - size(ROIframe,1)/2);
xCrop = round(xBest - size(ROIframe,2)/2);

ROIframesOut = frames((yCrop+10):(yCrop + size(ROIframe,1)-10), (xCrop+10):end, :);
