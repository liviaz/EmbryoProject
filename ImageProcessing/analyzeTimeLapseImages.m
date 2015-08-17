% analyze EEVA images of embryos 
% have outcome data of blast formation, and aneuploidy
% 
% Livia Zarnescu
% 3-25-15
%

% load outcome params
load('C:\Users\Livia\Desktop\Human embryo map\OutcomeParams.mat');
videoPath = 'C:\Users\Livia\Desktop\Human embryo map\Embryo movies\Processed';
addpath(videoPath);
addpath('C:\Users\Livia\Desktop\Human embryo map\');

%% 1. Extract features from each video time series

% Go through all embryos with outcome data, load in video, and extract
% params. First initialize param cell array

embryoParams = {};

for i = 1:length(dishNum)
    
    currFileName = ['BE' num2str(dishNum(i)) '_' wellNum{i} '.avi'];
    embryoParams{i} = extractEmbryoParams(videoPath, currFileName);
    
end

save('C:\Users\Livia\Desktop\Human embryo map\EmbryoParams.mat', ...
    'embryoParams');


%% 2. Save to .mat file

% lots of features per embryo




%% 2. Feed into SVM

% choose which features to use











