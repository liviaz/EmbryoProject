function out = getNorthPosition(desiredWidth, desiredHeight)
% out = getNorthPosition(desiredWidth, desiredHeight)
% Uses screensize to calculate distances from left and bottom
% out = [leftSpace bottomSpace width height]
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 

set(0,'Units','pixels');
screensize = get(0,'ScreenSize');
screenWidth = screensize(3);
screenHeight = screensize(4);
if screenWidth > desiredWidth
    leftSpace = floor((screenWidth - desiredWidth)/2);
    width = desiredWidth;
else
    leftSpace = 0;
    width = screenWidth;
end

extraSpace = 100;
if screenHeight > desiredHeight + extraSpace
    bottomSpace = screenHeight - (desiredHeight + extraSpace);
    height = desiredHeight;
else
    bottomSpace = 0;
    height =screenHeight-extraSpace;
end

out = [leftSpace bottomSpace width height];