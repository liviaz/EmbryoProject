function [] = runSURF(Acrop, ROIframe, fig1)

% zoom in on image area
figure(fig1);
h = helpdlg('Select ROI with pipette only');
pause(.5);
close(h);
coord = getrect;
coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
Apip = Acrop(coords(3):coords(4), coords(1):coords(2));

p1 = detectSURFFeatures(Apip, 'MetricThreshold', 20);
p2 = detectSURFFeatures(ROIframe, 'MetricThreshold', 20);

[feat1, vpts1] = extractFeatures(Apip, p1);
[feat2, vpts2] = extractFeatures(ROIframe, p2);

indexPairs = matchFeatures(feat1, feat2, 'MatchThreshold', 90);
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));

figure;
showMatchedFeatures(Apip, ROIframe, matchedPoints1, matchedPoints2);
legend('Matched Points 1', 'Matched Points 2');



