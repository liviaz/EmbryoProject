% analyze oocyte mechanics vs fertilization
% Livia Zarnescu Yanez
% 4-12-16


params = loadDataOocytes();
experimentType = params.experimentType;
numExperiments = params.numExperiments;
dateList = params.dateList;
numOocytes = params.numOocytes;
oocyteNums = params.oocyteNums;
fertInfo = params.fertInfo;
blastForm = params.blastForm;
hatchInfo = params.hatchInfo;
maturationEnv = params.maturationEnv;
k1ScaleFactor = params.k1ScaleFactor;
morphologyInfo = params.morphologyInfo;
measHour = params.measHour;

procDataPath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse Oocyte\';

% init param vectors
fertList = [];
blastList = [];
hatchList = [];
k0list = [];
k1list = [];
n0list = [];
n1list = [];
taulist = [];
colorMat = [];
aList = [];
matEnv = [];
legendParam = [];
oocyteM = [];
measHourList = [];
currE = 0;


exptsToPlot = [1 1 1 1 1 1 1 1 1 1 1 1];
aText = [];
cText = {};

% load in params by participant and embryo num
for i = 1:numExperiments
    if exptsToPlot(i) && experimentType(i)
        
        currDate = dateList{i}
        currDateU = strrep(currDate,'-','_');
        
        for j = 1:numOocytes(i)
            
            currE = currE + 1;
            currDataPath = [procDataPath currDate ' analysis\AutoMeasure\aspiration_data_' ...
                currDateU '_E' num2str(oocyteNums{i}(j)) '_16.mat'];
            
            if ~exist(currDataPath, 'file')
                currDataPath = [procDataPath currDate ' analysis\aspiration_data_' ...
                    currDateU '_E' num2str(oocyteNums{i}(j)) '_16.mat'];
            end
            
            if ~exist(currDataPath, 'file')
                currDataPath = [procDataPath currDate ' analysis\aspiration_data_' ...
                    currDateU '_E' num2str(oocyteNums{i}(j)) '.mat'];
            end
            
            if ~exist(currDataPath, 'file')
                currDataPath = [procDataPath currDate ' analysis\AutoMeasure\aspiration_data_' ...
                    currDateU '_E' num2str(oocyteNums{i}(j)) '.mat'];
            end
            
            % save embryo params and color
            if exist(currDataPath, 'file') && (morphologyInfo{i}(j) > -1) && (maturationEnv{i}(j) == 0)
                
                load(currDataPath);
                
                fertList = [fertList fertInfo{i}(j)];
                blastList = [blastList blastForm{i}(j)];
                hatchList = [hatchList hatchInfo{i}(j)];
                k0list = [k0list k0];
                k1list = [k1list (k1+k1ScaleFactor(i))];
                n0list = [n0list n0];
                n1list = [n1list n1];
                taulist = [taulist tau];
                oocyteM = [oocyteM morphologyInfo{i}(j)];
                matEnv = [matEnv maturationEnv{i}(j)];
                aText = [aText currE];
                measHourList = [measHourList measHour(i)];
                cText = {cText{:}, [currDate, '\_E', num2str(oocyteNums{i}(j))]};
                
                aPad = padarray(A,[0,65-length(A)], 'replicate', 'post');
                aList = [aList; aPad(1:65)];
                
                % color code based on fertilization and blast formation
                legendParamNums = 3;
                if (fertList(end) == 0)
                    legendParam = [legendParam 1];
                    colorMat = [colorMat; [.85 .65 .2]]; % no fert
                elseif (hatchList(end) == 0)
                    legendParam = [legendParam 2];
                    colorMat = [colorMat; [.9 .2 .8]]; % fert, no blast
                else
                    legendParam = [legendParam 3];
                    colorMat = [colorMat; [0 .2 .8]]; % fert, blast
                end
                

%                 % colorcode based on fertilization time
%                 legendParamNums = 3;
%                 if measHour(i) == 12
%                     legendParam = [legendParam 1];
%                     colorMat = [colorMat; [0 .2 .8]];%[.85 .65 .2]]; % 12 hrs
%                 elseif measHour(i) == 16
%                     legendParam = [legendParam 2];
%                     colorMat = [colorMat; [.5 .5 .5]]; % 16 hrs
%                 elseif measHour(i) == 20
%                     legendParam = [legendParam 3];
%                     colorMat = [colorMat; [.8 .2 .8]]; % 20 hrs
%                 end 
                    

                % colorcode based on maturation environment
%                 legendParamNums = 3;
%                 if (maturationEnv{i}(j) == 0)
%                     legendParam = [legendParam 1];
%                     colorMat = [colorMat; [.2 .6 .9]]; % in vivo
%                 elseif (maturationEnv{i}(j) == 1)
%                     legendParam = [legendParam 3];
%                     colorMat = [colorMat; [.85 .65 .2]]; % KSOM
%                 else
%                     legendParam = [legendParam 2];
%                     colorMat = [colorMat; [.9 .2 .8]]; % MM
%                 end
                
                % colorcode based on oocyte morphology
%                 legendParamNums = 3;
%                 if (morphologyInfo{i}(j) == 2)
%                     legendParam = [legendParam 1];
%                     colorMat = [colorMat; [.2 .6 .9]]; % MII
%                 elseif (morphologyInfo{i}(j) == 1)
%                     legendParam = [legendParam 2];
%                     colorMat = [colorMat; [.9 .2 .8]]; % MI
%                 else
%                     legendParam = [legendParam 3];
%                     colorMat = [colorMat; [.85 .65 .2]]; % GV
%                 end
                
            else
                colorMat = [colorMat; [NaN NaN NaN]];
                fertList = [fertList NaN];
                blastList = [blastList NaN];
                hatchList = [hatchList NaN];
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                n1list = [n1list NaN];
                taulist = [taulist NaN];
                oocyteM = [oocyteM NaN];
                matEnv = [matEnv NaN];
                legendParam = [legendParam NaN];
                aList = [aList; NaN*zeros(1,65)];
                aText = [aText NaN];
                measHourList = [measHourList NaN];
                cText = {cText{:}, ''};
            end
        end
    end
end

aList = 10^6*aList;
p1 = k1list;
p2 = n1list;

figure(1);
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
        xPatch = cos(theta)*(max(p1) - min(p1))/50 + p1(i);
        yPatch = sin(theta)*(max(p2) - min(p2))/50 + p2(i);

        h{i} = patch('xdata', xPatch, 'ydata', yPatch, 'edgecolor', colorMat(i,:), ...
            'facecolor', colorMat(i,:), 'facealpha', .6);

        hold on;
        if isempty(hLegend{legendParam(i)})
            i
            hLegend{legendParam(i)} = h{i};
        end
    end
end

set(gca, 'FontSize', 14);
title('Oocyte mechanics vs development');
grid on;
axis([.06 .145 -0.5 11]);
xlabel('k_1 parameter');
ylabel('\eta_1 parameter');
legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
    {'No fert', 'Fert', 'Hatched blast'}, 'Location', 'NorthEast');
% legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
%     {'12 hr', '16 hr', '20 hr'}, 'Location', 'NorthEast');
% legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
%     {'In vivo', 'MM', 'KSOM'}, 'Location', 'NorthEast');
% legend([hLegend{1}, hLegend{2}, hLegend{3}], ...
%     {'MII', 'MI', 'GV'}, 'Location', 'NorthEast');

b = num2str(aText');
c = cellstr(b);
dx = -0.01; dy = 0.015; dz = .5; % displacement so the text does not overlay the data points
hold on;
% text(p1+dx, p2+dy, cText'); %, p3+dz, c);
% view(.5,90)
% view(.5,0)


%% Plot effect of maturation environment on mechanics

figure(2);
clf;

% need to add .0314 to k1 for every .1psi increase in meas pressure
% need to add .4059 to n1 for every .1psi increase in meas pressure

k10 = bar(1, mean(k1list(matEnv == 0)), .8, 'facecolor', [.2 .6 .9]);
hold on;
ek10 = terrorbar(1, mean(k1list(matEnv == 0)), std(k1list(matEnv == 0)), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

k11 = bar(2, mean(k1list(matEnv == 2)), .8, 'facecolor', [.9 .2 .8]);
ek11 = terrorbar(2, mean(k1list(matEnv == 2)),std(k1list(matEnv == 2)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k12 = bar(3, mean(k1list(matEnv == 1)), .8, 'facecolor', [.85 .65 .2]);
ek12 = terrorbar(3, mean(k1list(matEnv == 1)),std(k1list(matEnv == 1)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);



set(gca, 'xtick', [1 2 3])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'in vivo', 'MM', 'KSOM'});
ylabel('k_1 parameter');
title('Maturation environment affects oocyte mechanics');
xlim([0.5 3.5]);

[p h] = ranksum(k1list(matEnv == 0), k1list(matEnv == 1))
[p h] = ranksum(k1list(matEnv == 1), k1list(matEnv == 2))


%% Plot correlation between morphology and mechanical params

figure(3);
clf;

k10 = bar(1, mean(k1list(oocyteM == 0)), .8, 'facecolor', [.2 .6 .6]);
hold on;
ek10 = terrorbar(1, mean(k1list(oocyteM == 0)), std(k1list(oocyteM == 0)), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

k11 = bar(2, mean(k1list(oocyteM == 1)), .8, 'facecolor', [.4 .8 .8]);
ek11 = terrorbar(2, mean(k1list(oocyteM == 1)),std(k1list(oocyteM == 1)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k12 = bar(3, mean(k1list(oocyteM == 2)), .8, 'facecolor', [.8 1 1]);
ek12 = terrorbar(3, mean(k1list(oocyteM == 2)),std(k1list(oocyteM == 2)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2 3])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'GV', 'MI', 'MII'});
ylabel('k_1 parameter');
title('Effect of nuclear maturation on mechanics');
xlim([0.5 3.5]);

[p h] = ranksum(k1list(oocyteM == 2), k1list(oocyteM == 1))
[p h] = ranksum(k1list(oocyteM == 1), k1list(oocyteM == 0))


%% Plot bar chart with all combinations of maturation environment and morphology


figure(4);
clf;

% In Vivo matured
k00 = bar(0, mean(k1list(oocyteM == 0 & matEnv == 0)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek00 = terrorbar(0, mean(k1list(oocyteM == 0 & matEnv == 0)), std(k1list(oocyteM == 0 & matEnv == 0)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

k10 = bar(.25, mean(k1list(oocyteM == 1 & matEnv == 0)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek10 = terrorbar(.25, mean(k1list(oocyteM == 1 & matEnv == 0)), std(k1list(oocyteM == 1 & matEnv == 0)), .1);
set(ek10, 'color', 'k', 'linewidth', 2);

k20 = bar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), .25, 'facecolor', [.8 1 1]);
hold on;
ek20 = terrorbar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), std(k1list(oocyteM == 2 & matEnv == 0)), .1);
set(ek20, 'color', 'k', 'linewidth', 2);

% MM matured
k02 = bar(1, mean(k1list(oocyteM == 0 & matEnv == 2)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek02 = terrorbar(1, mean(k1list(oocyteM == 0 & matEnv == 2)), std(k1list(oocyteM == 0 & matEnv == 2)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);

k12 = bar(1.25, mean(k1list(oocyteM == 1 & matEnv == 2)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek12 = terrorbar(1.25, mean(k1list(oocyteM == 1 & matEnv == 2)), std(k1list(oocyteM == 1 & matEnv == 2)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

k22 = bar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), .25, 'facecolor', [.8 1 1]);
hold on;
ek22 = terrorbar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), std(k1list(oocyteM == 2 & matEnv == 2)), .1);
set(ek22, 'color', 'k', 'linewidth', 2);

% KSOM matured
k01 = bar(2, mean(k1list(oocyteM == 0 & matEnv == 1)), .25, 'facecolor', [.2 .6 .6]);
hold on;
ek01 = terrorbar(2, mean(k1list(oocyteM == 0 & matEnv == 1)), std(k1list(oocyteM == 0 & matEnv == 1)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

k11 = bar(2.25, mean(k1list(oocyteM == 1 & matEnv == 1)), .25, 'facecolor', [.4 .8 .8]);
hold on;
ek11 = terrorbar(2.25, mean(k1list(oocyteM == 1 & matEnv == 1)), std(k1list(oocyteM == 1 & matEnv == 1)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k21 = bar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), .25, 'facecolor', [.8 1 1]);
hold on;
ek21 = terrorbar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), std(k1list(oocyteM == 2 & matEnv == 1)), .1);
set(ek21, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [0.25 1.25 2.25])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'in vivo', 'MM', 'KSOM'});
ylabel('k_1 parameter');
title('Effect of nuclear maturation on mechanics');
xlim([-0.5 3]);
%

legend([k01 k11 k21], {'GV', 'MI', 'MII'});

11111111
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0), k1list(oocyteM == 2 & matEnv == 1))
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 0), k1list(oocyteM == 2 & matEnv == 2))
[p h] = ranksum(k1list(oocyteM == 2 & matEnv == 1), k1list(oocyteM == 2 & matEnv == 2))




%% Bar plot of MII oocytes matured in vivo, MM, or KSOM


figure(5);
clf;

% In Vivo matured
k20 = bar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), .8, 'facecolor', [.8 1 1]);
hold on;
ek20 = terrorbar(.5, mean(k1list(oocyteM == 2 & matEnv == 0)), std(k1list(oocyteM == 2 & matEnv == 0)), .1);
set(ek20, 'color', 'k', 'linewidth', 2);

% MM matured
k22 = bar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), .8, 'facecolor', [.8 1 1]);
hold on;
ek22 = terrorbar(1.5, mean(k1list(oocyteM == 2 & matEnv == 2)), std(k1list(oocyteM == 2 & matEnv == 2)), .1);
set(ek22, 'color', 'k', 'linewidth', 2);

% KSOM matured
k21 = bar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), .8, 'facecolor', [.8 1 1]);
hold on;
ek21 = terrorbar(2.5, mean(k1list(oocyteM == 2 & matEnv == 1)), std(k1list(oocyteM == 2 & matEnv == 1)), .1);
set(ek21, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5])
set(gca, 'fontsize', 14);
ylim([0 .3]);
set(gca, 'xticklabel', {'in vivo', 'MM', 'KSOM'});
ylabel('k_1 parameter');
title('MII oocyte mechanics depend on environment');
xlim([0 3]);
%



%% Bar plot k1 of oocytes for non-fert/fert/blast


figure(5);
clf;
paramToPlot = n1list;


% Non-ferts
k00 = bar(.5, mean(paramToPlot(fertList == 0 & matEnv == 0)), .8, 'facecolor', [.85 .65 .2]);
hold on;
ek00 = terrorbar(.5, mean(paramToPlot(fertList == 0 & matEnv == 0)), std(paramToPlot(fertList == 0 & matEnv == 0)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% Fert, no blast
k01 = bar(1.5, mean(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0)), .8, 'facecolor', [.9 .2 .8]);
hold on;
ek01 = terrorbar(1.5, mean(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0)), ...
    std(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

% Fert + hatched blast
k02 = bar(2.5, mean(paramToPlot(hatchList == 1 & matEnv == 0)), .8, 'facecolor', [0 .2 .8]);
hold on;
ek02 = terrorbar(2.5, mean(paramToPlot(hatchList == 1 & matEnv == 0)), ...
    std(paramToPlot(hatchList == 1 & matEnv == 0)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5])
set(gca, 'fontsize', 14);
ylim([0 12]);
set(gca, 'xticklabel', {'non-fert', 'fert', 'hatch'});
ylabel('k_1 parameter');
title('MII oocyte mechanics correlated with development');
xlim([0 3]);
%

[p h] = ranksum(paramToPlot(fertList == 0 & matEnv == 0), paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0))
[p h] = ranksum(paramToPlot(fertList == 0 & matEnv == 0), paramToPlot(hatchList == 1 & matEnv == 0))
[p h] = ranksum(paramToPlot(fertList == 1 & hatchList == 0 & matEnv == 0), paramToPlot(hatchList == 1 & matEnv == 0))


%% Bar plot k1 of in vivo matured oocytes over time


figure(6);
clf;
paramToPlot = n1list;


% 12hr
k00 = bar(.5, mean(paramToPlot(measHourList == 12 & matEnv == 0)), .8, 'facecolor', [0 .2 .8]);
hold on;
ek00 = terrorbar(.5, mean(paramToPlot(measHourList == 12 & matEnv == 0)), std(paramToPlot(measHourList == 12 & matEnv == 0)), .1);
set(ek00, 'color', 'k', 'linewidth', 2);

% 16 hr
k01 = bar(1.5, mean(paramToPlot(measHourList == 16 & hatchList == 0 & matEnv == 0)), .8, 'facecolor', [.6 .6 .6]);
hold on;
ek01 = terrorbar(1.5, mean(paramToPlot(measHourList == 16 & hatchList == 0 & matEnv == 0)), ...
    std(paramToPlot(measHourList == 16 & hatchList == 0 & matEnv == 0)), .1);
set(ek01, 'color', 'k', 'linewidth', 2);

% 20 hr
k02 = bar(2.5, mean(paramToPlot(measHourList == 20 & matEnv == 0)), .8, 'facecolor', [.8 .2 .8]);
hold on;
ek02 = terrorbar(2.5, mean(paramToPlot(measHourList == 20 & matEnv == 0)), ...
    std(paramToPlot(measHourList == 20 & matEnv == 0)), .1);
set(ek02, 'color', 'k', 'linewidth', 2);


set(gca, 'xtick', [0.5 1.5 2.5])
set(gca, 'fontsize', 14);
ylim([0 12]);
set(gca, 'xticklabel', {'12 hr', '16 hr', '20 hr'});
ylabel('n_1 parameter');
title('MII oocyte mechanics depend on time');
xlim([0 3]);
%

[p h] = ranksum(paramToPlot(measHourList == 12 & matEnv == 0), paramToPlot(measHourList == 16 & matEnv == 0))
[p h] = ranksum(paramToPlot(measHourList == 16 & matEnv == 0), paramToPlot(measHourList == 20 & matEnv == 0))
[p h] = ranksum(paramToPlot(measHourList == 12 & matEnv == 0), paramToPlot(measHourList == 20 & matEnv == 0))



















