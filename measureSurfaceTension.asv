% measure surface tension of oocytes

date = '10-2-15';
fileDir = ['C:\Users\Livia\Desktop\IVF\Raw Data\Videos\Mouse Oocyte\videos ' date];
procDir = ['C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\' date ' analysis'];
type = 'Membrane';
eNum = 5;
embryoNum = ['E' num2str(eNum)];
convToMicrons = 40/108;
load([procDir '\' type '\all' type '.mat']);

pressureList = [];
depthList = [];

% fixed scale
pipCurrOpeningPixels = 63;
pipRefOpeningPixels = 126;


while (true)
    
    [fileName, pathName] = uigetfile([fileDir '\' embryoNum '*.jpg']);
    
    if fileName == 0
        break;
    end
    
    display(['Now measuring ' fileName]);
    
    A = rgb2gray(imread([pathName fileName]));
    f1 = figure(1);
    clf;
    imshow(A);
    movegui(f1, [-100, -100]);
    
    % record current pressure in psi
    pressure = input('Enter pressure (psi): ');

    % detect pipette edge
    % check if pipRef exists in fileDir
    if ~exist([fileDir '\pipRef.mat'], 'file')
        % save template of original
        MakeRoiTemplate(fileDir, .35);
    end
    
    % template matching from mat file
    load([fileDir '\pipRef.mat']);

    % find difference in scale (optional)
    figure(f1);
    h = helpdlg('Select ROI pipette only');
    pause(.5); close(h);
    coord = getrect;
    coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
    Apip = A(coords(3):coords(4), coords(1):coords(2));
    Acrop = A((coords(3)-20):(coords(4)+20), (coords(1)-20):(coords(2)+20));
    
    % to find scale conversion between screenshot and pixels
%     figure(f1);
%     clf;
%     imshow(Apip, 'InitialMagnification', 500);
%     h = helpdlg('Select inner corners of pipette');
%     pause(.5); close(h);
%     [x,y] = ginput(2);
%     
%     pipCurrOpeningPixels = abs(round(y(2) - y(1)));
    imScale = pipCurrOpeningPixels / pipRefOpeningPixels;
    ROIframeScaled = imresize(ROIframe, imScale);

    % now do template matching on scaled version 
    [~,I_NCC,~] = template_matching(ROIframeScaled, Acrop);
    [yBest,xBest] = find(I_NCC == max(I_NCC(:)));
    yCrop = round(yBest - size(ROIframeScaled,1)/2);
    xCrop = round(xBest - size(ROIframeScaled,2)/2);
    AcropPip = Acrop((yCrop+10):(yCrop + size(ROIframeScaled,1)-10), (xCrop+10):end, :);
    
    figure(f1);
    clf;
    imshow(AcropPip, 'InitialMagnification', 500);
    movegui(f1, [-100,-100]);
    
    % finally, select the depth into the pipette for the current pressure
    h = helpdlg('Select aspiration depth into pipette');
    pause(.5); close(h);
    [x, y] = ginput(1);
    
    if x < 0
        x = 0;
    end

    % save pressure and depth (convert back to pixels first)
    pressureList = [pressureList pressure];
    depthList = [depthList x/imScale*convToMicrons];
    
end

depthAll{eNum} = depthList;
pressureAll{eNum} = pressureList;

% plot
% close all;
figure(5);
clf;
scatter(pressureList, depthList, 150, 'marker', '+', 'linewidth', 3, ...
    'markeredgecolor', colorList(eNum,:));
set(gca, 'fontsize', 14);
xlabel('pressure (psi)');
ylabel('depth (\mum)');
title('Embryo mechanical properties');

% % polyfit
% p = polyfit(pressureList, depthList, 1);
% depthFit = polyval(p, pressureList);
% 
% hold on;
% plot(pressureList, depthFit, 'linewidth', 2, 'color', 'k');

% depthList is in MICRONS
save([procDir '\' type '\' embryoNum '_' lower(type) '.mat'], ...
    'pressureList', 'depthList');
save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll');



%%

% load settings
type = 'Zona';
load([procDir '\' type '\all' type '.mat']);

f6 = figure(6);
clf;
set(gca, 'fontsize', 14);
xlabel('stress (psi)');
ylabel('aspiration depth (\mum)');
title('zona measurement');
colorList = [0 0 .6; 0 .6 0; .6 0 0; 0 .6 .6; .6 0 .6];

% cellSize = zeros(1, length(pressureAll));
% surfTension = zeros(1, length(pressureAll));
% eqPressure = zeros(1, length(pressureAll));

plotHandles = [];

for i = 1:length(pressureAll)
    
    if ~isempty(pressureAll{i})
        
        hold on;
        subplot(2,3,i);
        h = scatter(pressureAll{i} - min(pressureAll{i}), ...
            depthAll{i}, 150, 'marker', '+', ...
            'linewidth', 3, 'markeredgecolor', colorList(i,:));
        set(gca, 'fontsize', 14);
        xlabel('stress (psi)');
        ylabel('aspiration depth (\mum)');
        title(['E' num2str(i)]);
        grid on;
        xlim([-.02, max(pressureAll{i}) + .02]);
        ylim([0 40]);

        plotHandles = [plotHandles h];
        
%         p = polyfit(pressureAll{i}, depthAll{i}, 1);
%         depthFit = polyval(p, pressureAll{i});
        
%         hold on;
%         plot(pressureAll{i}, depthFit, 'linewidth', 2, 'color', 'k');
        
%         [fileName, pathName] = uigetfile([fileDir '\E' num2str(i) '*.jpg']);
%         A = rgb2gray(imread([pathName fileName]));
%         figure(f6);
%         clf;
%         imshow(A);
%         movegui(f6, [-100, -100]);
%     
%         % find cell radius if measuring membrane / zona
%         diamList = zeros(1,2);
%         
%         for j = 1:2
%             [x,y] = ginput(2);
%             d = sqrt((x(1) - x(2))^2 + (y(1) - y(2))^2);
%             diamList(j) = d/imScale*convToMicrons;
%         end
%         
%         cellSize(i) = mean(diamList);

%         % find pressure for Lp = Rp, and subtract "holding" pressure
%         % eqPressure is in psi, cellSize is cell radius in microns when
%         % cell is aspirated at equipibrium (Lp = Rp)
%         pInv = polyfit(depthAll{i}, pressureAll{i}, 1);
%         pEq = polyval(pInv, 20) - min(pressureAll{i}); 
%         eqPressure(i) = pEq;
%         
%         % units of Pa * m = N / m
%         surfTension(i) = pEq*6895/(2*(1/(20*10^-6) - 1/(cellSize(i)*10^-6)));

    end
end

% xlim([-.05, .42]);
% ylim([0 .27]);
grid on;
% legend(plotHandles, {'E1', 'E2', 'E3', 'E4', 'E5'});

if isequal(type, 'Membrane')
    save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll', ...
    'cellSize', 'surfTension', 'eqPressure');
else
    save([procDir '\' type '\all' type '.mat'], 'depthAll', 'pressureAll', ...
    'cellSize');
end




