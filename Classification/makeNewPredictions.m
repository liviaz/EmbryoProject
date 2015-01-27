%
%   Inputs: testIndList tells which embryos are to be used in the training
%                       set to build the classifier (testIndList = 1) 
%                       and for which ones the morphology needs to be
%                       predicted (testIndList = 2)
%
%           mOut        has the ground truth data for all the embryos to be
%                       used in the training set and an arbitrary value
%                       (like 5 or something) for embryos to be predicted
%       
%           paramsOut   is the list of parameters for all the embryos (all
%                       the embryos in the training and test sets together)
%
%           numParams   is the number of parameters in the classifier
%
%           nGoodSelect tells the function to return the "nGoodSelect" best
%                       embryos as predicted by the classifier
%
%           nBadSelect  tells the function to return the "nBadSelect" worst
%                       embryos
%
%   Outputs: worstEmbryos   are the indices of the nBadSelect worst embryos
%                           as predicted by the classifier
%
%            bestEmbryos    are the indices of the nGoodSelect best embryos 
%                           as predicted by the classifier
%   
%            m_predict      are the predicted viability values of all the
%                           newly measured embryos
%
%
%            decDist        is the distance from the decision boundary for
%                           each embryo in the test set, gives a measure of
%                           the confidence of the predicted viability for
%                           each
%            decDistTrain   is the distance from the decision boundary for
%                           each embryo in the training set in order to
%                           predict likelihood of surviving for new embryos

function [worstEmbryos, bestEmbryos, m_predict, decDist, decDistTrain] = ...
    makeNewPredictions(testIndList, mOut, paramsOut, ...
    nGoodSelect, nBadSelect, fig_handle, enumList)

train = (testIndList == 1);
test = ~train;

mTrain = mOut(train);
mTest = mOut(test);

% the text below is already built into paramsOut

% % define list of parameters for embryos in each set
% % each is an numEmbryos x numParams matrix
% if numParams == 4
%     paramsTrain = paramsOut(train,:);
%     paramsTest = paramsOut(test,:);
% elseif numParams == 3
%     paramsTrain = paramsOut(train,1:3);
%     paramsTest = paramsOut(test,1:3);
% elseif numParams == 2
%     paramsTrain = paramsOut(train,1:2);
%     paramsTest = paramsOut(test,1:2);
% elseif numParams == 1
%     paramsTrain = paramsOut(train,1);
%     paramsTest = paramsOut(test,1);
% end

% optimize SVM params
posGroup = paramsOut(mOut == 1,:)';
negGroup = paramsOut(mOut == 0,:)';

zIn = OptimizeSvmParams(posGroup, ...
    negGroup, ...
    ones(1,size(posGroup,2)), ...
    ones(1,size(negGroup,2)), ...
    10, 'rbf', log([1 .5]), [5 5]);

paramsTrain = paramsOut(train,:);
paramsTest = paramsOut(test,:);

% make classifier and find predicted viability
%     figure(fig_handle_2);
embryoClassifier = svmtrain(paramsTrain, mTrain, 'kernel_function', ...
    'rbf', 'rbf_sigma', zIn(2), 'boxconstraint', zIn(1), 'method', 'smo');

% calculate predicted viability for test set
[m_predict, decDist] = svmclassify(embryoClassifier, paramsTest, ...
    'showplot', false);

% calculate distances from decision boundary for embryos in training set
[~, decDistTrain] = svmclassify(embryoClassifier, paramsTrain);

% highlight embryos furthest from decision boundary
[~, inds] = sort(decDist, 'ascend');
worstEmbryos = inds(end-nBadSelect+1:end)';
bestEmbryos = inds(1:nGoodSelect)';

% plot them along with embryos in training set
if size(paramsOut,2) > 2
    
    % plot all embryos in training set
    figure(fig_handle);
    ColorMat = zeros(length(mTrain), 3);
    num_pos = length(mTrain(mTrain == 1));
    ColorMat(mTrain == 1, :) = repmat([0 .6 0], num_pos, 1);
    ColorMat(mTrain == 0, :) = repmat([0 0 .6], length(mTrain) - num_pos, 1);
    hold on;
    scatter3(paramsTrain(:,1), paramsTrain(:,2), paramsTrain(:,3), ...
        100, ColorMat, 'filled');
    
    % overlay test embryos
    figure(fig_handle);
    hold on;
    % predicted worst
    scatter3(paramsTest(inds(end-nBadSelect+1:end),1), ...
        paramsTest(inds(end-nBadSelect+1:end),2), ...
        paramsTest(inds(end-nBadSelect+1:end),3), 100, [0 .6 .6], 'filled');
    hold on;
    % predicted best
    scatter3(paramsTest(inds(1:nGoodSelect),1), ...
        paramsTest(inds(1:nGoodSelect),2), ...
        paramsTest(inds(1:nGoodSelect),3), 100, [.7 .7 0], 'filled');
    hold on;
    % middle
    scatter3(paramsTest(inds(nGoodSelect+1:end-nBadSelect), 1), ...
        paramsTest(inds(nGoodSelect+1:end-nBadSelect), 2), ...
        paramsTest(inds(nGoodSelect+1:end-nBadSelect), 3), 100, [0 0 .6]);
    grid on;
    view(152, 20);
    
    hold on;
        
    a = enumList';
    b = num2str(a)
    c = cellstr(b);
    dx = .05; dy = 0.05; dz = .05; % displacement so the text does not overlay the data points
    text(paramsTest(a,1)+dx, paramsTest(a,2)+dy, paramsTest(a,3)+dz, c);
    
    
elseif size(paramsOut,2) == 2
    
    % plot all embryos in training set
    figure(fig_handle);
    ColorMat = zeros(length(mTrain), 3);
    num_pos = length(mTrain(mTrain == 1));
    ColorMat(mTrain == 1, :) = repmat([0 .6 0], num_pos, 1);
    ColorMat(mTrain == 0, :) = repmat([0 0 .6], length(mTrain) - num_pos, 1);
    hold on;
    scatter(paramsTrain(:,1), paramsTrain(:,2), ...
        100, ColorMat, 'filled');
    
    % overlay test embryos
    figure(fig_handle);
    hold on;
    % predicted worst
    h = scatter(paramsTest(inds(end-nBadSelect+1:end),1), ...
        paramsTest(inds(end-nBadSelect+1:end),2), ...
        100, [0 .6 .6], 'filled');
    hAxis = get(h,'parent');
    hold on;
    % predicted best
    scatter(paramsTest(inds(1:nGoodSelect),1), ...
        paramsTest(inds(1:nGoodSelect),2), ...
        100, [.7 .7 0], 'filled');
    hold on;
    % middle
    scatter(paramsTest(inds(nGoodSelect+1:end-nBadSelect), 1), ...
        paramsTest(inds(nGoodSelect+1:end-nBadSelect), 2), ...
        100, [0 0 .6]);
    grid on;
    
    hold on;
    a = enumList';
    b = num2str(a);
    c = cellstr(b);
    dx = .01; dy = 0.01; % displacement so the text does not overlay the data points
    text(paramsTest(a,1)+dx, paramsTest(a,2)+dy, c);
    
    hold on;
    plotSVandDC(embryoClassifier, hAxis, 1, 2, 0);
    
end