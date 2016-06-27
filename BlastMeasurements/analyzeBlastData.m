% analyze blast data

%% 1. load in data

blastList = [1 2 3 4 5];
repList = [3 4 3 5 4];
dataDir = 'C:\Users\Livia\Desktop\IVF\Raw Data\BlastMeasurements\MeasuredData';

numList = [];
k1list = [];
n1list = [];
colorMat = [];
currNum = 0;
legendParam = [];
aText = [];
cText = {};


for i = 1:length(blastList)
    for j = 1:repList(i)
   
        currNum = currNum + 1;
        currFileName = ['B' num2str(blastList(i)) '_loc' num2str(j)];
        
        if exist([dataDir '\' currFileName '.mat'], 'file')
            load([dataDir '\' currFileName '.mat']);
        
            k1list = [k1list k1];
            n1list = [n1list n1];
            numList = [numList i];
            aText = [aText currE];
            cText = {cText{:}, ['B' num2str(i) '\_L', num2str(j)]};

            
            legendParamNums = 5;
            if blastList(i) == 1
                legendParam = [legendParam 1];
                colorMat = [colorMat; [.2 .6 .8]];
            elseif blastList(i) == 2
                legendParam = [legendParam 2];
                colorMat = [colorMat; [.9 .2 .8]];
            elseif blastList(i) == 3
                legendParam = [legendParam 3];
                colorMat = [colorMat; [.85 .65 .2]];
            elseif blastList(i) == 4
                legendParam = [legendParam 4];
                colorMat = [colorMat; [.5 .5 .5]];
            elseif blastList(i) == 5
                legendParam = [legendParam 5];
                colorMat = [colorMat; [.2 .7 .2]];
            end
            
        end
            
    end
end

%% 2. bar plot

figure(2);
clf;
paramToPlot = n1list;

% Blast 1
k00 = bar(.5, mean(paramToPlot(numList == 1)), .8, 'facecolor', [.5 .5 .5]);
hold on;
ek00 = terrorbar(.5, mean(paramToPlot(numList == 1)), std(paramToPlot(numList == 1)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% Blast 2
k01 = bar(1.5, mean(paramToPlot(numList == 2)), .8, 'facecolor', [.5 .5 .5]);
hold on;
ek01 = terrorbar(1.5, mean(paramToPlot(numList == 2)), ...
    std(paramToPlot(numList == 2)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

% Blast 3
k02 = bar(2.5, mean(paramToPlot(numList == 3)), .8, 'facecolor', [.5 .5 .5]);
hold on;
ek02 = terrorbar(2.5, mean(paramToPlot(numList == 3)), ...
    std(paramToPlot(numList == 3)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);

% Blast 4
k03 = bar(3.5, mean(paramToPlot(numList == 4)), .8, 'facecolor', [.5 .5 .5]);
hold on;
ek03 = terrorbar(3.5, mean(paramToPlot(numList == 4)), ...
    std(paramToPlot(numList == 4)), .1);
set(ek03, 'color', 'k', 'linewidth', 2);

% Blast 5
k03 = bar(4.5, mean(paramToPlot(numList == 5)), .8, 'facecolor', [.5 .5 .5]);
hold on;
ek03 = terrorbar(4.5, mean(paramToPlot(numList == 5)), ...
    std(paramToPlot(numList == 5)), .1);
set(ek03, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5 3.5 4.5])
set(gca, 'fontsize', 14);
% ylim([0 12]);
set(gca, 'xticklabel', {'B1', 'B2', 'B3', 'B4', 'B5'});
ylabel('n_1 parameter');
title('trophectoderm mechanics');
xlim([0 5]);
%


%% 3. scatter plot

p1 = k1list;
p2 = n1list;


figure(3);
clf;
hold on;
h = cell(1,length(k1list));
hLegend = cell(1,legendParamNums);

for i = 1:length(k1list)
    i
    if ~isnan(legendParam(i))
        i
%         h{i} = plot(p1(i), p2(i), 'marker', 'o', 'markeredgecolor', ...
%             colorMat(i,:), 'markersize', 12, 'color', 'none', ...
%             'linewidth', 4);

        theta = 0:0.1:(2*pi);     
        xPatch = cos(theta)*(max(p1) - min(p1))/30 + p1(i);
        yPatch = sin(theta)*(max(p2) - min(p2))/25 + p2(i);

        h{i} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorMat(i,:), ...
            'facecolor', colorMat(i,:), 'facealpha', .7);

        hold on;
        if isempty(hLegend{legendParam(i)})
            i
            hLegend{legendParam(i)} = h{i};
        end
    end
end

set(gca, 'FontSize', 14);
title('trophectoderm mechanics');
grid on;
% axis([.06 .145 -0.5 11]);
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');
legend([hLegend{1}, hLegend{2}, hLegend{3}, hLegend{4}, hLegend{5}], ...
    {'B1', 'B2', 'B3', 'B4', 'B5'}, 'Location', 'NorthEast');


b = num2str(aText');
c = cellstr(b);
dx = -0.01; dy = 0.015; % displacement so the text does not overlay the data points
hold on;
% text(p1+dx, p2+dy, cText'); 
