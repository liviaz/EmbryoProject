% load all data for running classification using SVM
% Inputs:
%   inputMethod: if inputMethod = 1, output training and test set by
%                                    randomly splitting data evenly
%
%                if inputMethod = 2, output training set as all data so
%                                    far, and test data as current
%                                    experiment (choose this option for
%                                    real-time classification)
%
%                                    the embryos to be classified must have
%                                    morphology values of 5, while all the
%                                    rest must have either 2 or 4
%
%   nGroups: Number of partitions to split data for validation
%            If inputMethod = 2, nGroups MUST equal 2
%            If inputMethod = 1, nGroups must be >= 2
%
% Outputs:
%   mOut = ground truth morphology data (contains test + train data)
%   paramsOut = param lists for test + train data
%   testIndList = indices indicating which measurement is in which group
%                 If inputMethod = 2, then testIndList = 1 is the training
%                       data and testIndList = 2 is the testing data
%                 If inputMethod = 1, then each value of testIndList
%                       represents a different partition of the data


function [mOut, paramsOut, testIndList, enumList] = loadSVMdataOocyte(inputMethod, ...
    nGroups, plotInput, filePath1, paramNumsToUse)

if plotInput
    figPlot = figure;
end

% ========================================
% MAKE SURE LIST OF EXPERIMENTS IS CORRECT
% ========================================


% whichToPlot(1) is 5-14-13 (fresh)
% whichToPlot(2) is 5-30-13 (frozen)
% whichToPlot(3) is 2-11-14 (fresh)
% whichToPlot(4) is 2-27-14 (fresh)
whichToPlot = [1 1 1 1];
plotAuto = [ones(1,length(whichToPlot))];
pixelConv = [108*ones(1,length(whichToPlot))];

% ===== MAKE SURE NEW MORPHOLOGY DATA IS SAVED =====
% =========== FILL IN HERE =========================
saveNewMorphology('mouse oocyte');
load('morphologyMouseOocyte.mat');

% morphology results of all experiments so far
% =========== FILL IN HERE ======================
morphology = cell(1,length(whichToPlot));
morphology{1} = morphology5_14_13;
morphology{2} = morphology5_30_13;
morphology{3} = morphology2_11_14;
morphology{4} = morphology2_27_14;

% time lapse parameters of all experiments so far
% =========== FILL IN HERE ======================
timeLapse = cell(1,length(whichToPlot));
timeLapse{1} = timeLapse5_14_13;
timeLapse{2} = timeLapse5_30_13;
timeLapse{3} = timeLapse2_11_14;
timeLapse{4} = timeLapse2_27_14;

% dates of all experiments so far
% =========== FILL IN HERE ======================
dateList = cell(1,length(whichToPlot));
dateList{1} = '5-14-13';
dateList{2} = '5-30-13';
dateList{3} = '2-11-14';
dateList{4} = '2-27-14';

dateList2 = cell(1,length(whichToPlot));
for j = 1:length(whichToPlot)
    dateU = dateList{j};
    dateUI = strfind(dateU, '-');
    dateU(dateUI) = '_';
    dateList2{j} = dateU;
end

% for making legend to scatter plot
colorMat = [];
k0list = [];
k1list = [];
n0list = [];
taulist = [];
elonglist = [];
n1list = [];
morphologyOut = [];
timeLapseOut = [];

for jNum = 1:length(whichToPlot)
    if whichToPlot(jNum)
        
        cCurr = zeros(length(morphology{jNum}), 3);
        s = size(cCurr(morphology{jNum} == 2, :));
        cCurr(morphology{jNum} == 2, :) = repmat([0 0 .6], s(1), 1);
        s = size(cCurr(morphology{jNum} == 3, :));
        cCurr(morphology{jNum} == 3, :) = repmat([0 .6 .6], s(1), 1);
        s = size(cCurr(morphology{jNum} == 4, :));
        cCurr(morphology{jNum} == 4, :) = repmat([0 .6 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 1, :));
        cCurr(morphology{jNum} == 1, :) = repmat([1 0 0], s(1), 1);
        s = size(cCurr(morphology{jNum} == 5, :));
        cCurr(morphology{jNum} == 5, :) = repmat([.7 .7 0], s(1), 1);
        
        colorMat = [colorMat ; cCurr];
        
        % add on to morphologyOut
        % 2 = didn't fertilize
        % 4 = fertilized
        morphologyAdd = morphology{jNum};
        morphologyAdd(morphologyAdd > 1) = 4;
        morphologyAdd(morphologyAdd < 2) = 2;

        % add to morphologyOut
        morphologyOut = [morphologyOut morphologyAdd];        
        
        % add to timeLapseOut
        if ~isempty(timeLapse{jNum})
            timeLapseOut = [timeLapseOut; timeLapse{jNum}];
        else
            timeLapseOut = [timeLapseOut; NaN*ones(length(morphology{jNum}),3)];
        end        
        
        for iNum = 1:length(morphology{jNum})
            
            embryoString = num2str(iNum);
            
            % if plotAuto is on, see if an automatically measured file exists.
            % If not, just use the manually measured one.
            % if plotAuto is off, just use the manually measured one
            testPath = [filePath1, dateList{jNum}, ' analysis\AutoMeasure\aspiration_data_', ...
                dateList2{jNum}, '_E', embryoString, '.mat'];
            
            if plotAuto(jNum) && exist(testPath)
                filePath2 = [dateList{jNum} ' analysis\AutoMeasure\aspiration_data_' ...
                    dateList2{jNum}];
            else
                filePath2 = [dateList{jNum} ' analysis\aspiration_data_' ...
                    dateList2{jNum}];
            end
            
            if exist([filePath1 filePath2 '_E', embryoString, '.mat']) && ...
                    ~isnan(morphology{jNum}(iNum))
                
                load([filePath1 filePath2 '_E', embryoString, '.mat']);
                
                aspiration_depth = aspiration_depth * 40 * 10^-6 / pixelConv(jNum); % convert from pixels to meters
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                taulist = [taulist tau];
                n1list = [n1list n1];
                elonglist = [elonglist F0/(k0 + k1)];
                
                currColor = cCurr(iNum,:);
                xdata = t;
                ydata = aspiration_depth;
                
                if plotInput
                    figure(figPlot);
                    plot(xdata, ydata, 'ob', 'Color', currColor);
                    hold on;
                    plot(xfine - min(xfine), yfit, 'Color', currColor);
                end
                
            else
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                taulist = [taulist NaN];
                n1list = [n1list NaN];
                elonglist = [elonglist NaN];
            end
        end
        
    end
end

if plotInput
    xlim([0 .45]);
    set(gca, 'FontSize', 14);
    xlabel('time (seconds)');
    ylabel('aspiration depth (\mum)');
    title('Aspiration Depth of Mouse Oocytes');
end

% put all variables and labels together in a list
% get rid of NaNs (morphologyOut has some NaNs, and so does k1list, etc
% but in different locations

% if not using any time lapse params, can use all mechanical data
if isempty(intersect(paramNumsToUse,[5 6 7]))
    numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
        ~isnan(morphologyOut);
else
    numsToPlot = ~isnan(k1list) & ~isnan(n1list) & ~isnan(taulist) & ...
        ~isnan(morphologyOut) & ~isnan(timeLapseOut(:,1))';
end

if inputMethod == 2
    numToTest = length(morphologyOut(morphologyOut == 5));
    enumList = 1:numToTest;
    enumList = enumList(numsToPlot(end-numToTest+1:end));
else
    enumList = [];
end

% new names for variables
k = k1list(numsToPlot);
n = n1list(numsToPlot);
t = taulist(numsToPlot);
k2 = k0list(numsToPlot);

% n = log(n);
% k2 = log(k2);

% normalize each parameter list from 0-1
kk = (k - min(k)) / (max(k) - min(k));
nn = (n - min(n)) / (max(n) - min(n));
tt = (t - min(t)) / (max(t) - min(t));
kk2 = (k2 - min(k2)) / (max(k2) - min(k2));

% Scale data first
% normalize all params 0-1
m = morphologyOut(numsToPlot);
m(m == 2) = 0;
m(m == 4) = 1;

% [k1, n1, tau, k0, TL1, TL2, TL3]
paramsOut = [kk ; nn ; tt ; kk2]';

% process time lapse parameters and add to paramsOut
if ~isempty(intersect(paramNumsToUse,[5 6 7]))
    
    timeLapseOut = timeLapseOut(numsToPlot,:);
    
    timeLapseOut(:,1) = (timeLapseOut(:,1) - min(timeLapseOut(:,1))) / ...
        (max(timeLapseOut(:,1)) - min(timeLapseOut(:,1)));
    timeLapseOut(:,2) = (timeLapseOut(:,2) - min(timeLapseOut(:,2))) / ...
        (max(timeLapseOut(:,2)) - min(timeLapseOut(:,2)));
    timeLapseOut(:,3) = (timeLapseOut(:,3) - min(timeLapseOut(:,3))) / ...
        (max(timeLapseOut(:,3)) - min(timeLapseOut(:,3)));
    
    paramsOut = [paramsOut timeLapseOut];
end

paramsOut = paramsOut(:,paramNumsToUse);
mOut = m;

% ============================================================size(para
% BEGINNING OF IF STATEMENT
% ============================================================

if inputMethod == 1
    
    % Create a n-fold cross-validation to compute classification error
    % m is ground truth data (+/-1)
    testIndList = crossvalind('Kfold', m, nGroups)';
    
elseif inputMethod == 2
    
    testIndList = ones(1,length(mOut));
    testIndList(mOut > 1) = 2;
    
end

end


