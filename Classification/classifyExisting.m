%
% if crossValidate = 1, predictOut and decDistOut are the same length as
% mOut. If crossValidate = 0, predictOut and decDistOut are only as long as
% the test group, as indicated by where testIndList = 2
%
%

function [predictOut, decDistOut, zOut] = classifyExisting(testIndList, mOut, ...
    paramsOut, fig_handle, crossValidate, nGroups, plotInput, optimFlag, zIn)

% ======================================================
% Perform cross-validation between the n groups
% ======================================================

% close all;
% create class performance object
cpEmbryos = classperf(mOut, 'positive', 1, 'negative', 0);

% sigma is 1
% C is .001
% fig_handle_2 = figure;

if crossValidate
    
    % vector to save all predictions in
    predict_all = zeros(1,length(mOut));
    decDist_all = zeros(1,length(mOut));
    
    for i = 1:nGroups
        
        if optimFlag && i == 1
            
            % just set a value for rbf_sigma and C, don't optimize
            
                        posGroup = paramsOut(mOut == 1,:)';
                        negGroup = paramsOut(mOut == 0,:)';
            
%                         zOut = OptimizeSvmParams(posGroup, ...
%                             negGroup, ...
%                             ones(1,size(posGroup,2)), ...
%                             ones(1,size(negGroup,2)), ...
%                             nGroups, 'rbf', log([1 .6]), [5 5]);
                        
%                         zOut = OptimizeSvmParams(posGroup, ...
%                             negGroup, ...
%                             ones(1,size(posGroup,2)), ...
%                             ones(1,size(negGroup,2)), ...
%                             nGroups, 'rbf', [1 .6], [2 2]);
            
            zOut = [1 .7];
            zIn = zOut;
            rbf_sigma = zOut(2);
            C = zOut(1);
            
        else
            rbf_sigma = zIn(2);
            C = zIn(1);
            zOut = zIn;
        end

        test = (testIndList == i);
        train = ~test;
        
        mTrain = mOut(train);
        mTest = mOut(test);
        
        % define list of parameters for embryos in each set
        % each is an nEmbryos x numParams matrix
        paramsTrain = paramsOut(train, :);
        paramsTest = paramsOut(test, :);
        
        % make classifier and find predicted viability
        %     figure(fig_handle_2);
        
        embryoClassifier = svmtrain(paramsTrain, mTrain, 'kernel_function', ...
            'rbf', 'rbf_sigma', rbf_sigma, 'boxconstraint', C, ...
            'showplot', false, 'method', 'smo', 'autoscale', true);
        [m_predict, decDist] = svmclassify(embryoClassifier, paramsTest, ...
            'showplot', false);
        predict_all(test) = decDist < 0;
        decDist_all(test) = decDist;
       
        
        if plotInput
            
            if size(paramsOut,2) < 3
                axisLims = [0 1 0 1 0 1];
            else
                axisLims = [min(paramsOut(:,1)) max(paramsOut(:,1)) ...
                    min(paramsOut(:,2)) max(paramsOut(:,2)) ...
                    min(paramsOut(:,3)) max(paramsOut(:,3))];
            end
            
            % plot results
            plotClassificationResults(mTest, decDist, paramsTest, fig_handle, ...
                embryoClassifier, axisLims);
        end
        
    end
    
%     % account for NaNs
%     mOut
%     predict_all
%     mOut(isnan(predict_all)) = NaN;
    
    true_pos = length(mOut(mOut == 1 & predict_all == 1));
    false_pos = length(mOut(mOut == 0 & predict_all == 1));
    true_neg = length(mOut(mOut == 0 & predict_all == 0));
    false_neg = length(mOut(mOut == 1 & predict_all == 0));
    
    sensitivity = true_pos / (true_pos + false_neg)
    specificity = true_neg / (true_neg + false_pos)
    ppv = true_pos / length(predict_all(predict_all == 1))
    
%     classperf(cpEmbryos, predict_all) 
    predictOut = predict_all;
    decDistOut = decDist_all;
    
else
    
    % if not doing cross-validation, just do the classification once 
    test = (testIndList == 1);
    train = ~test;
    
    mTrain = mOut(train);
    mTest = mOut(test);
    
    C = zIn(1);
    rbf_sigma = zIn(2);
    zOut = zIn;
    
    % define list of parameters for embryos in each set
    % each is an nx3 matrix
    paramsTrain = paramsOut(train, :);
    paramsTest = paramsOut(test, :);
    
    % make classifier and find predicted viability
    %     figure(fig_handle_2);
    embryoClassifier = svmtrain(paramsTrain, mTrain, 'kernel_function', ...
        'rbf', 'rbf_sigma', rbf_sigma, 'boxconstraint', C, ...
        'showplot', false);
    [m_predict, decDist] = svmclassify(embryoClassifier, paramsTest, ...
        'showplot', false);
    
    if plotInput
        % plot results
        plotClassificationResults(mTest, decDist, paramsTest, fig_handle, ...
            embryoClassifier);
    end
    
    size(mTest)
    size(m_predict)
    
    true_pos = length(mTest(mTest == 1 & m_predict' == 1));
    false_pos = length(mTest(mTest == 0 & m_predict' == 1));
    true_neg = length(mTest(mTest == 0 & m_predict' == 0));
    false_neg = length(mTest(mTest == 1 & m_predict' == 0));
    
    sensitivity = true_pos / (true_pos + false_neg)
    specificity = true_neg / (true_neg + false_pos)
    ppv = true_pos / length(m_predict(m_predict == 1))
    
    classperf(cpEmbryos, m_predict', test) 
    predictOut = m_predict';
    decDistOut = decDist;
    
end