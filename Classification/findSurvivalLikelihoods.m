% given predicted survival values and distances from classifier decision
% boundary, 

function survivalLikelihoods = findSurvivalLikelihoods(decDistTest, ...
    decDistTrain, mOut, testIndList, methodToUse)

train = (testIndList == 1);
test = ~train;

mTrain = mOut(train);
mTest = mOut(test);

% decDistTest corresponds to mTest
% decDistTrain corresponds to mTrain

if methodToUse == 1
    
    % ========================================================================
    % Method 1: Estimate based on how many viable and nonviable embryos have a
    % similar distance to the decision boundary
    % ========================================================================
    
    % for a given embryo in the test set, calculate its likelihood of survival
    % by looking at the relative numbers of viable/non-viable embryos in the
    % training set with decision boundary distances within +/- decInterval
    
    % want at least 10 embryos in each "bin"
    decInterval = (max([decDistTest' decDistTrain']) - ...
        min([decDistTest' decDistTrain'])) / (length(decDistTrain') / 10);
    
    survivalLikelihoods = zeros(1,length(mTest));
    
    for i = 1:length(mTest)
        
        currDist = decDistTest(i);
        
        % find embryos in training set with decision distances close to currDist
        closeDistEmbryos = mTrain(decDistTrain < currDist + decInterval & ...
            decDistTrain > currDist - decInterval);
        
        numPos = length(closeDistEmbryos(closeDistEmbryos == 1));
        numNeg = length(closeDistEmbryos(closeDistEmbryos == 0));
        
        survivalLikelihoods(i) = numPos / (numPos + numNeg);
        
    end

elseif methodToUse == 2
    
    % ========================================================================
    % Method 2: Estimate based on Platt's method, fit to logistic curve
    % ========================================================================    
    
        
%         minFun = inline(['sum(mTest.*(AB(1)*decDistTest + AB(2)) ' ...
%         '+ log(1 + exp(-AB(1)*decDistTest - AB(2))))'], ...
%         'mTest', 'decDistTest', 'AB');
    
    [ABout fval] = fminsearch(@(AB) findSigmoidParams(AB, mTrain', decDistTrain), [0 0]);
    
    survivalLikelihoods = ((1 + exp(ABout(1)*decDistTest + ABout(2))).^(-1))';
    
end
