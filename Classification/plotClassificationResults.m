% plot results of embryo classification
% mTest is the ground truth data in the test data set
% m_predict is the predicted embryo viability as output by the classifier
% paramsTest is a nx3 matrix with the parameters of the embryos in the test
% set
% fig_handle is just the handle to the figure to plot in
% inputMethod = 1 if ground truth is known and classifier is being
% cross-validated. InputMethod = 2 if making predictions for new embryos

function plotClassificationResults(mTest, decDist, paramsTest, fig_handle, ...
    embryoClassifier, axisLims)

m_predict = (decDist < -.5);

% Plot Results for Test data set
% green and blue for actual
ColorMat = zeros(length(mTest), 3);
num_pos = length(mTest(mTest == 1));
ColorMat(mTest == 1, :) = repmat([0 .6 0], num_pos, 1);
ColorMat(mTest == 0, :) = repmat([0 0 .6], length(mTest) - num_pos, 1);

% yellow and light blue for predicted ones
ColorPredict = zeros(length(m_predict), 3);
num_pos = length(m_predict(m_predict == 1));
ColorPredict(m_predict' == 1, :) = repmat([.7 .7 0], num_pos, 1);
ColorPredict(m_predict' == 0, :) = repmat([0 .6 .6], length(m_predict) - num_pos, 1);


if size(paramsTest,2) > 2
    
%     hold on;
%     subplot(1, 2, 1);
    % make 2 scatter plots
    figure(fig_handle);
    hold on;
    h = scatter3(paramsTest(:,1), paramsTest(:,2), paramsTest(:,3), ...
        100, ColorMat, 'filled');
    set(gca, 'FontSize', 14);
    set(h, 'Marker', 'o');
    title('Actual Viability');
%     xlabel('k1 parameter');
%     ylabel('n1 parameter');
%     zlabel('tau parameter');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    zlim([min(paramsTest(:,3)) max(paramsTest(:,3))+.1]);
    view(152,20);
    grid on;
    axis(axisLims);
    
    hC = flipud(get(h, 'children'));
    
    legend([hC(find(mTest == 1, 1, 'first')) ...
        hC(find(mTest == 0, 1, 'first'))], ...
        'Viable', 'Nonviable', 'Location', 'North');

%     hold on;
%     a = currEmbryoNums';
%     b = num2str(a);
%     c = cellstr(b);
%     dx = .05; dy = 0.05; dz = .05; % displacement so the text does not overlay the data points
%     text(paramsTest(:,1)+dx, paramsTest(:,2)+dy, paramsTest(:,3)+dz, c);
    
    figure(fig_handle+1);
    hold on;
%     subplot(1, 2, 2);
%     hold on;
    h = scatter3(paramsTest(:,1), paramsTest(:,2), paramsTest(:,3), ...
        100, ColorPredict, 'filled');   
    set(gca, 'FontSize', 14);
    set(h, 'Marker', 'o');
    title('SVM Predicted Viability');
%     xlabel('k1 parameter');
%     ylabel('n1 parameter');
%     zlabel('tau parameter');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    zlim([min(paramsTest(:,3)) max(paramsTest(:,3))+.1]);
    view(152,20);
    grid on;
    axis(axisLims);
    
    hC = flipud(get(h, 'children'));
    
    legend([hC(find(m_predict == 1, 1, 'first')) ...
        hC(find(m_predict == 0, 1, 'first'))], ...
        'Predicted Viable', 'Predicted NonViable', 'Location', 'North');
    
elseif size(paramsTest,2) == 2
    
    figure(fig_handle);
    hold on;
    h = scatter(paramsTest(:,1), paramsTest(:,2), ...
        100, ColorMat, 'filled');
    hAxis = get(h,'parent');
    set(h, 'Marker', 'o');
    set(gca, 'FontSize', 14);
    title('Actual Viability');
%     xlabel('k1 parameter');
%     ylabel('\tau parameter');
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    axis(axisLims(1:4));

    hold on;
    plotSVandDC(embryoClassifier, hAxis, 1, 2, -.2);
    %         hold on;
    %         a = currEmbryoNums';
    %         b = num2str(a);
    %         c = cellstr(b);
    %         dx = .025; dy = .025; % displacement so the text does not overlay the data points
    %         text(paramsTest(:,1)+dx, paramsTest(:,2)+dy, c);
    
    figure(fig_handle+1);
    hold on;
    h = scatter(paramsTest(:,1), paramsTest(:,2), ...
        100, ColorPredict, 'filled');
    set(h, 'Marker', 'o');
    set(gca, 'FontSize', 14);
    title('SVM Predicted Viability');
%     xlabel('k1 parameter');
%     ylabel('\tau parameter');
    axis(axisLims(1:4));
    xlim([min(paramsTest(:,1)) max(paramsTest(:,1))+.1]);
    ylim([min(paramsTest(:,2)) max(paramsTest(:,2))+.1]);
    
    hold on;
    if length(get(gca, 'children')) == 5
        plotSVandDC(embryoClassifier, hAxis, 1, 2, 0); 
    end
    
end
