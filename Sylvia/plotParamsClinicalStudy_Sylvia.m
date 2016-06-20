% plot params for all patients in clinical study
% Livia Zarnescu Yanez 11-10-15



%% 1. Load in a single data file
% and plot it


filePath = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E1.mat'
load(filePath)

x1 = k0
y1 = k1

filePath2 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E2.mat'
load(filePath2)
x2 = k0
y2 = k1

filePath3 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E3.mat'
load(filePath3)

x3 = k0
y3 = k1

filePath4 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E4.mat'
load(filePath4)

x4 = k0
y4 = k1

filePath5 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E5.mat'
load(filePath5)
x5 = k0
y5 = k1

filePath6 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E6.mat'
load(filePath6)
x6 = k0
y6 = k1

filePath7 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E7.mat'
load(filePath7)
x7 = k0
y7 = k1

filePath8 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E8.mat'
load(filePath8)
x8 = k0
y8 = k1

filePath9 = '\Users\Livia\Desktop\IVF\Processed Data\Human\MECH001\MECH001_E9.mat'
load(filePath9)
x9 = k0
y9 = k1

figure
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9, 'marker','*')




%%

params = loadDataClinicalStudy_Sylvia();


baseFilePath = '\Users\Livia\Desktop\IVF\Processed Data\Human\';
patientName = 'MECH00';
numPatients = 22;
numEmbryos = [9 0 14 13 10 7 0 11 13 6 9 0 21 13 0 11 7 0 6 0 6 9]; % finish filling this in];


% these are the variables which will store all the loaded data
% they are initialized as empty vectors
k1list = [];
k0list = [];
n0list = [];
n1list = [];
viability = [];
colorMat = [];
ICSI = [];
outcomeInfo =[];
patientAge = [];
patientBMI = [];



% this is the for loop which will do the same thing multiple times
% it takes a variable i which changes value each time the loop is executed
for i = 1:numPatients
    for j = 1:numEmbryos(i)
        
        % this statement will print the value of i and j each time the loop is
        % executed
        fprintf('%f\n', i);
        fprintf('%f\n', j);
        
        
        if (i < 10)
            currPatientName = ['MECH00' num2str(i)]
        else
            currPatientName = ['MECH0' num2str(i)]
        end
        
        
        
        % step 1: construct the file path for the current embryo
        currFilePath = [baseFilePath currPatientName '\' ...
            currPatientName '_E' num2str(j) '.mat'] % complete this
        
        % step 2: load the file at that path
        load(currFilePath);
        
        
        % step 3: add the loaded data to a vector
        k1list = [k1list k1]; % and so on
        k0list = [k0list k0];
        n0list = [n0list n0];
        n1list = [n1list n1];
        
        outcomeInfo = [outcomeInfo params.outcomeInfo{i}(j)];
        ICSI = [ICSI params.ICSI{i}(j)]
        patientAge = [patientAge params.patientAge{i}]
        
        
        if params.outcomeInfo{i}(j) == 1
            if params.ICSI {i}(j) == 1
                colorMat = [colorMat; [.9 .5 0]];
            else
                colorMat = [colorMat; [0 .6 0]];
            end
        else
            if params.ICSI {i}(j) == 1
                colorMat = [colorMat; [0 0 .8]];
            else
                colorMat = [colorMat; [1 0 0]];
            end
        end
        
        
        
    end
end

% only make scatterplot of data with ICSI == 1



% now make a scatterplot of all the data you loaded in
close all;
figure;
scatter3(k1list, outcomeInfo, ICSI, 100, colorMat, 'filled');









%% 1. Load in data

params = loadDataClinicalStudy();
numParticipants = params.numParticipants;
participantIDs = params.participantIDs;
numEmbryos = params.numEmbryos;
outcomeInfo = params.outcomeInfo;
zygoteM = params.zygoteM;
d3M = params.d3M;
blastM = params.blastM;
ICSI = params.ICSI;
PGD = params.PGD;
gender = params.gender;
patientAge = params.patientAge;

totalNumEmbryos = sum(numEmbryos);
procDataPath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human\';

% init param vectors
mList = [];
k0list = [];
k1list = [];
n0list = [];
n1list = [];
taulist = [];
gList = [];
colorMat = [];
icsiList = [];
zList = [];
aList = [];
ptList = [];
currE = 0;
dayThreeM1 = [];
dayThreeM2 = [];
blastStage = [];
blastScore = [];

ptsToPlot = [ones(1,21) 1]; %[1 0 0 0 0 0 0 0 0 0 0 0 0];
aText = [];
cText = {};

% load in params by participant and embryo num
for i = 1:numParticipants
    if ptsToPlot(i)
        i
        for j = 1:numEmbryos(i)
            
            currE = currE + 1;
            currDataPath = [procDataPath participantIDs{i} '\' ...
                participantIDs{i} '_E' num2str(j) '.mat'];
            ptList = [ptList i];
            
            % save embryo params and color
            if exist(currDataPath, 'file') && (ICSI{i}(j)== 0)
                load(currDataPath);
                mList = [mList outcomeInfo{i}(j)>0];
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                n1list = [n1list n1];
                taulist = [taulist tau];
                gList = [gList gender{i}(j)];
                icsiList = [icsiList ICSI{i}(j)];
                zList = [zList zygoteM{i}(j)];
                aText = [aText currE];
                cText = {cText{:}, ['P', num2str(i), '\_E', num2str(j)]};
                %                 cText = {cText{:}, num2str(j)};
                
                aPad = padarray(A,[0,65-length(A)],'replicate', 'post');
                aList = [aList; aPad(1:65)];
                
                % get day 3 morphology data
                dayThreeM = d3M{i}{j};
                if numel(dayThreeM) == 0
                    dayThreeM1 = [dayThreeM1 0];
                    dayThreeM2 = [dayThreeM2 0];
                else
                    dayThreeM1 = [dayThreeM1 str2num(dayThreeM(1))];
                    dayThreeM2 = [dayThreeM2 getNumFromNumeral(dayThreeM(2:end))];
                end
                
                % get day 5 morphology data
                blastString = blastM{i}{j};
                if numel(blastString) == 0
                    blastStage = [blastStage 0];
                    blastScore = [blastScore 0];
                else
                    blastStage = [blastStage str2num(blastString(1))];
                    if numel(blastString) < 3
                        blastScore = [blastScore 0];
                    else
                        blastScore = [blastScore (double((blastString(2)-65)+(blastString(3)-65))/2)];
                    end
                end
                
                
                if (blastStage(end) > 4) && (blastScore(end) < 1.6)
                    colorMat = [colorMat; [.2 .6 .9]];
                elseif (blastStage(end) > 0)
                    colorMat = [colorMat; [.85 .65 .2]];%[.9 .2 .8]];%
                else
                    colorMat = [colorMat; [.85 .65 .2]];
                end
                
                %                 if outcomeInfo{i}(j) == 2 % unknown
                %                     mList(end) = 2;
                %                     colorMat(end,:) = [1 0 0];
                %                 end
                
                % implantation/pregnancy color coding
                if outcomeInfo{i}(j) == 10 % no hCG rise
                    mList(end) = 10;
                    colorMat(end,:) = [1 0 0];
                end
                
                if outcomeInfo{i}(j) == 11 % hCG rise, no pregnancy
                    mList(end) = 11;
                    colorMat(end,:) = [1 0 0];%[.7 .7 0];
                end
                
                if outcomeInfo{i}(j) == 12 % hCG rise, pregnancy
                    mList(end) = 12;
                    colorMat(end,:) = [0 .5 0];
                end
                
                % ICSI / IVF color coding
                %                     if ICSI{i}(j) == 1
                %                         colorMat(end,:) = [.2 .6 .9];
                %                     elseif ICSI{i}(j) == 0
                %                         colorMat(end,:) = [.85 .65 .2];
                %                     end
                
                % PGD color coding
                %                   if PGD{i}(j) == 1
                %                       colorMat(end,:) = [.2 .6 .9]; % euploid
                %                   elseif PGD{i}(j) == 0
                %                       colorMat(end,:) = [.85 .65 .2]; % aneuploid
                %                   else
                %                       k1list(end) = NaN;
                %                   end
                
                
                
            else
                colorMat = [colorMat; [NaN NaN NaN]];
                mList = [mList NaN];
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                n1list = [n1list NaN];
                taulist = [taulist NaN];
                icsiList = [icsiList NaN];
                aList = [aList; NaN*zeros(1,65)];
                gList = [gList NaN];
                zList = [zList NaN];
                aText = [aText NaN];
                cText = {cText{:}, ''};
                dayThreeM1 = [dayThreeM1 NaN];
                dayThreeM2 = [dayThreeM2 NaN];
                blastStage = [blastStage NaN];
                blastScore = [blastScore NaN];
            end
        end
    end
end


aList = 10^6*aList;
% colorMat(icsiList == 1,:) = repmat([NaN NaN NaN], length(icsiList(icsiList == 1)), 1);

transfer1 = k1list(mList==10);
transfer2 = n1list(mList==10);
transferC = colorMat(mList==10,:);

%% 2. Plot

% now plot!
p1 = taulist;
p2 = blastScore; %n1list;
% p3 = ptList; %zList + .1*randn(1,length(zList));


figure(1); clf;
% h = scatter3(p1, p2, p3, 200, colorMat, 'filled');
h = scatter(p1,p2,200,colorMat,'linewidth',4);
% h = scatter(p1(ptList<13 | mList==0), p2(ptList<13 | mList==0), 200, colorMat(ptList<13 | mList==0,:), 'linewidth', 4);
hold on;
%h2 = scatter(p1(ptList==13 & mList>0),p2(ptList==13 & mList>0),200,[1 0 0],'linewidth', 3);
%h3 = scatter(transfer1,transfer2,200,transferC,'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('Blastocyst morphology vs mechanics');
xlabel('tau parameter');
ylabel('blastScore');
% zlabel('zygote morphology');
% axis([min(p1) max(p1) min(p2) max(p2) min(p3) max(p3)]);
% set(gca, 'yscale', 'log')
% axis([min(p1) max(p1) min(p3) max(p3)]);
% set(gca, 'zscale', 'linear');
% grid on;
% xlim([.3 .8]);
% ylim([0 .8]);

b = num2str(aText');
c = cellstr(b);
dx = -0.01; dy = 0.015; dz = .5; % displacement so the text does not overlay the data points
hold on;
% text(p1+dx, p2+dy, cText'); %, p3+dz, c);
view(.5,90)


%% 1. Plot ICSI vs IVF mech properties

figure(2);
clf;
k11 = bar(1, mean(k1list(icsiList == 0)), .8, 'facecolor', [0 .6 .6]);
hold on;
% ek11 = errorbar(1, mean(k1list(icsiList == 0)),std(k1list(icsiList == 0)),'color', 'k', 'linewidth', 2);
ek11 = terrorbar(1,mean(k1list(icsiList == 0)), std(k1list(icsiList == 0)), .1);
set(ek11, 'color', 'k', 'linewidth', 2);

k12 = bar(1.8, mean(k1list(icsiList == 1)), .8, 'facecolor',[.6 .9 .9] );
ek12 = terrorbar(1.8, mean(k1list(icsiList == 1)),std(k1list(icsiList == 1)), .1);
set(ek12, 'color', 'k', 'linewidth', 2);

n11 = bar(3.2, mean(n1list(icsiList == 0)), .8, 'facecolor', [1 .7 .2]);
en11 = terrorbar(3.2, mean(n1list(icsiList == 0)),std(n1list(icsiList == 0)), .1);
set(en11, 'color', 'k', 'linewidth', 2);

n12 = bar(4, mean(n1list(icsiList == 1)), .8, 'facecolor',[1 .9 .6]);
en12 = terrorbar(4, mean(n1list(icsiList == 1)),std(n1list(icsiList == 1)), .1);
set(en12, 'color', 'k', 'linewidth', 2);

set(gca, 'xtick', [1 2])
set(gca, 'fontsize', 14);
ylim([0 1]);
set(gca, 'xtick', [1.4 3.6]);
set(gca, 'xticklabel', {'k_1', '\eta_1'});
% xlabel('fertilization method');
ylabel('parameter value');
% title('fertilization method affects mechanics');
xlim([0.5 4.5]);

[p h] = ranksum(k1list(icsiList == 0), k1list(icsiList == 1))
[p h] = ranksum(n1list(icsiList == 0), n1list(icsiList == 1))

%% 2. Plot all aspiration depth curves

tVec = (0:1:64)/60;
figure(3);
clf;


for i = 1:length(mList)
    if ~isnan(mList(i))
        plot(tVec, 10^6*aList(i,:), 'color', colorMat(i,:), 'marker', 'o');
        hold on;
    end
end

axis([0 max(tVec) 0 45]);

%% 4. Extract and plot morphology



figure(1); clf;
h = scatter(dayThreeM1 + .75*rand(1,length(dayThreeM1)), ...
    dayThreeM2  + .75*rand(1,length(dayThreeM1)), 200, colorMat, 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('Mechanical parameters predict blastocyst formation');
xlabel('day3 morphology parameter 1');
ylabel('day3 morphology parameter 2');
% zlabel('zygote morphology');
% axis([min(p1) max(p1) min(p2) max(p2) min(p3) max(p3)]);
% set(gca, 'yscale', 'log')
% axis([min(p1) max(p1) min(p3) max(p3)]);
% set(gca, 'zscale', 'linear');
grid on;
% xlim([.3 .7]);
% ylim([0 .8]);

% b = num2str(aText');
% c = cellstr(b);
% dx = -0.01; dy = 0.015; dz = .5; % displacement so the text does not overlay the data points
% hold on;
% text(p1+dx, p2+dy,p3+dz, cText'); %, p3+dz, c);
% view(.5,90)


%% 5. Bar chart comparison of morphology vs mechanics

figure;
% day 3 morphology ROC
h1 = bar(.85, 0.756, .3, 'facecolor', [.8 .3 .3]);
hold on;
e1 = errorbar(.85,0.756,.02, 'color', 'k', 'linewidth', 2);
setErrorBar(e1,.85,.03);

% mechanics ROC
h2 = bar(1.15, .848, .3, 'facecolor', [.3 .7 .8]);
e2 = errorbar(1.15,.848,.011, 'color', 'k', 'linewidth', 2);
setErrorBar(e2,1.15,.03);

% day 3 morphology PR
h3 = bar(1.85, 0.599, .3, 'facecolor', [.8 .3 .3]);
e3 = errorbar(1.85,0.599,.045, 'color', 'k', 'linewidth', 2);
setErrorBar(e3,1.85,.03);

% mechanics PR
h4 = bar(2.15, .703, .3, 'facecolor', [.3 .7 .8]);
e4 = errorbar(2.15, .703, .016, 'color', 'k', 'linewidth', 2);
setErrorBar(e4,2.15,.03);

set(gca, 'xtick', [1 2])
set(gca, 'fontsize', 14);
set(gca, 'xticklabel', {'ROC', 'PR'});
ylabel('Area under curve');
title('Comparison of mechanics and day 3 morphology');
xlim([0.5 2.5]);
grid on;

%% 6. SVM figure plot


figure(2);
svmStruct = svmtrain([k1list; n1list]', mList > 0, 'kernel_function', ...
    'rbf', 'rbf_sigma', 1, 'boxconstraint', .5, 'ShowPlot', true);

c = get(gca, 'children');
figure(3);
clf;

h = scatter(p1,p2,200,colorMat,'linewidth',4);
hold on;
contour(c(1).XData, c(1).YData, c(1).ZData, -.3*[1 1], 'k', 'linewidth', 2);
set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('ICSI');
xlabel('k1 parameter');
ylabel('n1 parameter');
grid on;
xlim([.3 .8]);
ylim([0 .8]);

close(2);







