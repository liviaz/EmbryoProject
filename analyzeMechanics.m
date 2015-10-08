% Livia Zarnescu
% 10-7-2015
% Make some comparisons between zona, membrane, and bulk props

clear all;
close all;

%% 1. Compare bulk properties at different pressures

% 10-2-15 data
date = '10-2-15';
eNums = [4 5];
pressureVals = [1 2 3 4];
baseDir = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte';
dataDir = [baseDir '\' date ' analysis\AutoMeasure'];
colorVec = [0 0 .6; 0 .6 0; .6 0 0; 0 .6 .6];
markerList = {'+', 'o', '^', '*'};

k1M = zeros(length(eNums),length(pressureVals));
n1M = zeros(length(eNums),length(pressureVals));
k0M = zeros(length(eNums),length(pressureVals));
tauM = zeros(length(eNums),length(pressureVals));

f1 = figure(1);
clf;

% plot aspiration curves for different pressures
for i = 1:length(eNums)
    for j = 1:length(pressureVals)
        
        currFile = [dataDir, '\aspiration_data_', strrep(date, '-', '_'), ...
            '_E', num2str(eNums(i)), '_pt', num2str(pressureVals(j)), '.mat'];

        load(currFile);
        
        k1M(i,j) = k1
        n1M(i,j) = n1;
        k0M(i,j) = k0;
        tauM(i,j) = tau;
        
        figure(1);
        hold on;
        scatter(t, A*10^6, 150, 'marker', markerList{j}, 'linewidth', 2, ...
            'markeredgecolor', colorVec(i,:));
        plot(xfine-min(xfine), yfit*10^6, 'linewidth', 2, 'color', colorVec(i,:));
        
    end
end

% plot how parameters change for different pressures
f2 = figure(2);
clf;

subplot(2,2,1);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], k1M')
colormap([0 0 .6; 0 .6 0]);
title('k1');
xlim([0 .5]);

subplot(2,2,2);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], n1M')
colormap([0 0 .6; 0 .6 0]);
title('n1');
xlim([0 .5]);

subplot(2,2,3);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], k0M')
colormap([0 0 .6; 0 .6 0]);
title('k0');
xlim([0 .5]);

subplot(2,2,4);
set(gca, 'fontsize', 14);
bar([.1 .2 .3 .4], tauM')
colormap([0 0 .6; 0 .6 0]);
title('tau');
xlim([0 .5]);


%% 2. correlate surface tension with zona stiffness and bulk props



































