%
%   Procedure for measuring embryos and making predictions
%
%   1. Run depthDetectionAll to measure all embryos
%
%   2. Edit saveNewMorphology to add in new experiment with value "5" for
%      all new (unknown ground truth) embryos
%
%   3. Edit loadSVMdata
%
%   4. Run this function with inputMethod = 2 and nGroups = 2
%
%
% inputMethod = 1 does classification on existing ground truth data to
% evaluate goodness of classifier

% inputMethod = 2 does classification of new measurements with unknown
% ground truth


function [AUC_ROC, AUC_PR, zSVM, aROCstd, aPRstd, Ytotal, Ztotal, decDistTest] = ...
    classifyEmbryos(inputMethod, nGroups, method, paramNumsToUse, numRepeats, plotInput)

% close all;
% clear all;

if nargin < 4
    paramNumsToUse = 1:3;
end
if nargin < 5
    numRepeats = 10;
end
if nargin < 6
    plotInput = 0;
end

addpath('C:\Users\Livia\Desktop\IVF\Code\Matlab Code');
addpath('C:\Users\Livia\Desktop\SVM code\Bruce Code\MG_ImpactDetection\Livia');
addpath('C:\Users\Livia\Desktop\SVM code\Bruce Code\MG_ImpactDetection');

if isequal(method, 'human embryo')
    filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human analysis\';
elseif isequal(method, 'mouse oocyte')
    filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse oocyte analysis\';
else % mouse embryo
    filePath1 = 'C:\Users\Livia\Desktop\IVF\Processed Data\Mouse embryo analysis\';
end

% inputMethod = 1;
% nGroups = 5;

Xref = 0:.01:1;
zSVM = [0 0];
Atotal = [];
AtotalPR = [];
decDistTest = [];
Ytotal = [];
Ztotal = [];
AUC_ROC = [];
AUC_PR = [];
aROCstd = [];
aPRstd = [];

if inputMethod == 1
    figROC = figure;
elseif inputMethod == 2
    numRepeats = 1;
end

if plotInput
    fig_handle = figure;
else
    fig_handle = NaN;
end


for j = 1:numRepeats
    
    j
    %% first load data to be used for classification
    
    if isequal(method, 'human embryo')
        saveNewMorphology('human');
        [mOut, paramsOut, testIndList, enumList] = loadSVMdataHuman(inputMethod, nGroups, 0, ...
            filePath1, paramNumsToUse);
    elseif isequal(method, 'mouse oocyte')
        saveNewMorphology('mouse oocyte');
        [mOut, paramsOut, testIndList, enumList] = loadSVMdataOocyte(inputMethod, nGroups, 0, ...
            filePath1, paramNumsToUse);
    else % mouse embryo
        saveNewMorphology('mouse embryo');
        [mOut, paramsOut, testIndList, enumList] = loadSVMdata(inputMethod, nGroups, 0, ...
            filePath1, paramNumsToUse);
    end
    
    
    
    %% Make predictions about test group
    
    nGoodSelect = 6;
    nBadSelect = 6;
    
    if inputMethod == 2
        
        % make predictions about newly measured embryos
        [worstEmbryos, bestEmbryos, m_predict, decDistTest, decDistTrain] = ...
            makeNewPredictions(testIndList, mOut, paramsOut, ...
            nGoodSelect, nBadSelect, fig_handle, enumList);%, [.4 1.6]);
        
        decDistTest
        
        worstEmbryos = sort(worstEmbryos)
        bestEmbryos = sort(bestEmbryos)
        
        % find likelihood of each "test" embryo surviving
        survivalLikelihoods = findSurvivalLikelihoods(decDistTest, ...
            decDistTrain, mOut, testIndList, 2);
        
        % find best and worst embryos with this method
        [survivalSort, sInd] = sort(survivalLikelihoods, 'ascend');
        worstEmbryos2 = sInd(1:nBadSelect);
        bestEmbryos2 = sInd(end-nGoodSelect+1:end);
        
        % finally, display
        fprintf('\n Embryos Least Likely to Survive:');
        worstEmbryos2
        
        fprintf('\n Survival Likelihoods of Least Viable Embryos:');
        survivalLikelihoods(worstEmbryos2)
        
        fprintf('\n Embryos Most Likely to Survive:');
        bestEmbryos2
        
        fprintf('\n Survival Likelihoods of Most Viable Embryos:');
        survivalLikelihoods(bestEmbryos2)
        
        zSVM = [0 0];
        
        
    elseif inputMethod == 1
        
        crossValidate = 1;
        
        % do cross-validation on data set that already has ground truth
        [~, decDistTest, zSVM] = classifyExisting(testIndList, mOut, ...
            paramsOut, fig_handle, crossValidate, nGroups, plotInput, (j == 1), zSVM);
        
        %% Plot ROC curves
        
        % AUC = .85
        % AUC of precision recall curve is .93
        if crossValidate
            [X, Y, T, AUC] = perfcurve(mOut, -decDistTest, 1);
            prevalence = length(mOut(mOut == 1)) / length(mOut);
        else
            mTrue = mOut(testIndList == 1);
            [X, Y, T, AUC] = perfcurve(mTrue, -decDistTest, 1);
            prevalence = length(mTrue(mTrue == 1)) / length(mTrue);
        end
        
        % get Yi out with points at a standard set of locations (at Xref)
        % otherwise finding average ROC curve isn't really possible
        Yi = interpForRoc(Xref, X, Y);
        Atotal = [Atotal AUC];
        
        if j == 1
            %         Xtotal = Xi/numRepeats; % X is 1 - specificity
            Ytotal = Yi/numRepeats; % Y is sensitivity
%             Atotal = AUC;
        else
            %         Xtotal = Xtotal + Xi/numRepeats;
            Ytotal = Ytotal + Yi/numRepeats;
%             Atotal = Atotal + AUC;
        end
        
        Zi = Yi*prevalence ./ (Yi*prevalence + Xref*(1-prevalence)); % PPV
        Yi = Yi(~isnan(Zi));
        Zi = Zi(~isnan(Zi));
        
        AtotalPR = [AtotalPR trapz([0 Yi 1], [Zi(1) Zi Zi(end)])];
        
        figure(figROC);
        hold on;
        plot(X,Y);
        hold on;
        
    end
    
    
end

% calculate and plot average ROC curve if doing cross-validation
if inputMethod == 1
    
    if 1%plotInput
        figure;
        set(gca, 'FontSize', 14);
        hold on;
        plot([Xref 1], [Ytotal 1], 'Color', [0 0 1], 'LineWidth', 3);
        xlabel('1 - Specificity'); ylabel('Sensitivity')
        title('Average ROC curve');
        grid on;
        axis([0 1 0 1]);
    end
    
    Ztotal = Ytotal*prevalence ./ (Ytotal*prevalence + Xref*(1-prevalence)); % PPV
    
    if 1%plotInput
        figure, plot(Ytotal, Ztotal, 'Color', [0 0 1], 'LineWidth', 3);
        set(gca, 'FontSize', 14);
        xlabel('Recall (sensitivity)'); ylabel('Precision (PPV)');
        title('Average PR curve');
        grid on;
        axis([0 1 0 1]);
    end
    
%     Ytotal = Ytotal(~isnan(Ztotal));
%     Ztotal = Ztotal(~isnan(Ztotal));
    
    
    AUC_ROC = mean(Atotal);
    aROCstd = std(Atotal);
    AUC_PR = mean(AtotalPR);
    aPRstd = std(AtotalPR);
%     AUC_ROC = Atotal / numRepeats;
%     aROCstd = NaN;
%     AUC_PR = trapz([0 Ytotal 1], [Ztotal(1) Ztotal Ztotal(end)]);
%     aPRstd = NaN;
    
    fprintf(['\nArea Under ROC curve is ' num2str(AUC_ROC) ...
        ' with stdev ' num2str(aROCstd) '\n']);
    fprintf(['\nArea Under PR curve is ' num2str(AUC_PR) ...
        ' with stdev ' num2str(aPRstd) '\n']);
    
else
    
    AUC_ROC = NaN;
    AUC_PR = NaN;
    aROCstd = NaN;
    aPRstd = NaN;
    
end

end
