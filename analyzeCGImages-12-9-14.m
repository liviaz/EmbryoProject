% analyze CG images 
% correlate with mechanical params
% 12-9-14

clear all;
close all;

%% 1. For each embryo, save tiff stack

% 5, 31, 34 are missing
imageDirectory = 'C:\Users\Livia\Desktop\IVF\CG stain\ZonaHardening12-9-14\';

% embryoNums = [1 3 4 5 6 7 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ...
%     26 27 28 29 31 32 33 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53];

% if we ignore the ones with zona remaining, we are left with
embryoNums = [4 5 7 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ...
    26 27 28 29 32 33 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53];

dirNames = cell(1, length(embryoNums));

for i = 1:length(embryoNums)
   dirNames{i} = ['E' num2str(embryoNums(i)) '\']; 
   fileNames{i} = ['E' num2str(embryoNums(i))];
end

dirNames = {dirNames, 'control\'};
fileNames = {fileNames, 'control'};
tifNames = fileNames;


fileToSave = 'data/CG_imaging_12_9_14.mat';

for i = 1:length(dirNames)
    
    display(['Filename: ' tifNames(i)]);
    
    % get list of all tif images in directory and 
    % make empty matrix to hold them
    imageList = dir([imageDirectory dirNames{i}]);
    imageList = imageList(3:end);
    imageProps = imfinfo([imageDirectory dirNames{i} imageList(1).name]);
    imageMatrix = zeros(imageProps.Width, imageProps.Height, length(imageList));
    
    for j = 1:length(imageList)
        currFileBaseName = [imageDirectory dirNames{i} imageList(j).name];
        A = imread(currFileBaseName);
        A = im2double(A(:,:,2));
        imageMatrix(:,:,j) = A;
    end
    
    figure, imagesc(mean(imageMatrix,3));
    colorbar;
    caxis([0 .15]);
    pause;
    
    WriteTiffStack(imageMatrix, [imageDirectory 'TIFS\' tifNames{i} '.tif']);
    
end


%% 2. Make struct containing ROIs for each embryo

% build struct
% sub-structs for group, and embryo within each group
groupNames = tifNames;
embryosInImage = ones(1,length(embryoNums));
embryosInGroup = ones(1,length(embryoNums));
imageGroup = 1:length(embryoNums);
eNumOffset = zeros(1,length(embryoNums)); % for groups with multiple TIFS
embryoInfo = struct;

for i = 1:length(groupNames)
    embryoInfo.(groupNames{i}) = struct;
    i
    for j = 1:sum(embryosInImage(imageGroup == i))
        embryoInfo.(groupNames{i}).(['E' num2str(j)]) = struct;
        j
        % this step needs to be fixed if using multiple embryos/image
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).tifName = tifNames{i};
    end
end

embryoInfo.embryosInImage = embryosInImage;
embryoInfo.groupNames = groupNames;
embryoInfo.imageGroup = imageGroup;
embryoInfo.tifNames = tifNames;
embryoInfo.imageDirectory = imageDirectory;
embryoInfo.embryosInGroup = embryosInGroup;
embryoInfo.eNumOffset = eNumOffset;
embryoInfo.embryoNums = embryoNums;

save(fileToSave, 'embryoInfo');

%% Read in Tiffs and extract ROIs

load(fileToSave);
embryosInImage = embryoInfo.embryosInImage;
groupNames = embryoInfo.groupNames;
eNumOffset = embryoInfo.eNumOffset;

for i = 1:length(embryoInfo.tifNames)
    
    currFileName = [embryoInfo.imageDirectory '\TIFS\' embryoInfo.tifNames{i} '.tif'];
    
    % calculate projection image
    imInfo = imfinfo(currFileName); % get dimensions of first image
    imageMatrix = zeros(imInfo(1).Height, imInfo(1).Width); % green channel
    
    for k = 1:length(imInfo)
        currImage = im2double(imread(currFileName, 'Index', k));
        imageMatrix = imageMatrix + currImage;
    end
    
    % take average and apply blur
    imageMatrix = imageMatrix / length(imInfo);
    imageMatrix = imfilter(imageMatrix, fspecial('gaussian', 10, 2));
    
    % display image
    figure(1);
    clf;
    imagesc(imageMatrix);
    set(gca, 'fontsize', 14);
    
    imageGroup = embryoInfo.imageGroup(i);
    
    % select ROI and save
    for j = 1:embryosInImage(i)
        
        figure(1);
        title(['draw region around embryo ' num2str(j)]);
        
        % second image will have second embryo from first group
        eNum = j + eNumOffset(i);
        display(['Tiff # ' num2str(i), ', embryo #' num2str(eNum)]);
            
        coords = getrect;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).coords = coords;
        ROI = imageMatrix(round(coords(2):(coords(2)+coords(4))), ...
            round(coords(1):(coords(1)+coords(3))));
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).ROI = ROI;
        
        figure(2);
        clf;
        imshow(imadjust(ROI));
        
    end
end

save(fileToSave, 'embryoInfo');

%% Detect cell boundary and cytoplasm

close all;
load(fileToSave);
embryosInImage = embryoInfo.embryosInImage;
groupNames = embryoInfo.groupNames;
eNumOffset = embryoInfo.eNumOffset;

for i = 1:length(embryoInfo.tifNames)
    
    imageGroup = embryoInfo.imageGroup(i);
    
    % select ROI and save
    for j = 1:embryosInImage(i)

        figure(1);
        title(['draw region around embryo ' num2str(j)]);
        
        % second image will have second embryo from first group
        eNum = j + eNumOffset(i);
        display(['Tiff # ' num2str(i), ', Group # ' ...
            num2str(imageGroup) ', embryo #' num2str(eNum)]);
            
        ROI = embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).ROI;
        
        figure(2);
        clf;
        imshow(imadjust(ROI));
        
        maxRadius = mean(size(ROI))/2;
        radiusRange = round(linspace(maxRadius*.7, maxRadius, 20));
        [mask, mR, mX, mY] = getCellMask(ROI, radiusRange);    
        
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mask = mask;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mR = mR;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mX = mX;
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).mY = mY;
                
    end
end

save(fileToSave, 'embryoInfo');


%% Extract cell params

close all;
load(fileToSave);
embryosInGroup = embryoInfo.embryosInGroup;
groupNames = embryoInfo.groupNames;

for i = 1:length(embryosInGroup)
    
    % select ROI and save
    for j = 1:embryosInGroup(i)
        
        ROI = embryoInfo.(groupNames{i}).(['E' num2str(j)]).ROI;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).group = i;
        
        figure(1);
        clf;
        imshow(imadjust(ROI));
        
         % get mean brightness of representative ROI
        coords = getrect;
        subROI = ROI(round(coords(2):(coords(2)+coords(4))), ...
            round(coords(1):(coords(1)+coords(3))));
        embryoInfo.(groupNames{imageGroup}).(['E' num2str(eNum)]).subROI = subROI;
        
        % now calculate params!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % param for overall cytoplasm brightness, maybe brightness going
        % radially outward, and maybe brightness around outside of cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [brightParam, circParam, radParam] = getROIparams_CGimages(embryoInfo, i, j);

        embryoInfo.(groupNames{i}).(['E' num2str(j)]).brightParam = brightParam;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam = circParam;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).radParam = radParam;
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).circStd = std(circParam);
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).circMean = mean(circParam);
        embryoInfo.(groupNames{i}).(['E' num2str(j)]).subROImean = mean(mean(subROI));
        
    end
end

save(fileToSave, 'embryoInfo');


%% Make figures from params

load(fileToSave);
embryosInGroup = embryoInfo.embryosInGroup;
groupNames = embryoInfo.groupNames;
circMeanList = cell(1,40);
circStdList = cell(1,40);
radMeanList = cell(1,40);
brightParamList = cell(1,40);
circDiffList = cell(1,40);
ROImeanList = cell(1,40);
groupList = cell(1,40);
figure(1);
clf;
hold on;

for i = 1:length(embryosInGroup)
    
    % accumulate param vectors
    for j = 1:embryosInGroup(i)
        
        groupList{i} = [groupList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).group];
        circMeanList{i} = [circMeanList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).circMean];
        circStdList{i} = [circStdList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).circStd];
        radMeanList{i} = [radMeanList{i} mean(embryoInfo.(groupNames{i}).(['E' num2str(j)]).radParam)];
        brightParamList{i} = [brightParamList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).brightParam];
        circDiffList{i} = [circDiffList{i} (max(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam) ...
            - min(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam))/...
            embryoInfo.(groupNames{i}).(['E' num2str(j)]).circMean];
        ROImeanList{i} = [ROImeanList{i} embryoInfo.(groupNames{i}).(['E' num2str(j)]).subROImean];
        
        if mod(i, 3) == 2
            currColor = [0 .6 0];
            h1 = plot(linspace(0, 2*pi, 100), ...
                smooth(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam), ...
                'color', currColor, 'linewidth', 2);
        else 
            currColor = [0 0 .6];
            h2 = plot(linspace(0, 2*pi, 100), ...
                smooth(embryoInfo.(groupNames{i}).(['E' num2str(j)]).circParam), ...
                'color', currColor, 'linewidth', 2);
        end
        
        
    
    end
end

set(gca, 'fontsize', 14);
xlabel('angle (rad)');
ylabel('image intensity (a.u.)');
title('Image Intensity Around Cell Boundary');
xlim([0 2*pi]);
legend([h1 h2], 'viable', 'nonviable');

embryoInfo.groupList = groupList;
embryoInfo.circMeanList = circMeanList;
embryoInfo.circStdList = circStdList;
embryoInfo.radMeanList = radMeanList;
embryoInfo.brightParamList = brightParamList;
embryoInfo.circDiffList = circDiffList;
embryoInfo.ROImeanList = ROImeanList;
save(fileToSave, 'embryoInfo');


%% plot CG data from 12-9-14
% modeled off of "analyzeMicroinjectionData.m"

fileToSave = 'data/CG_imaging_12_9_14.mat';
load(fileToSave);
embryoNums = [4 7 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ...
    26 27 28 29 32 33 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53];

mechNums = [4 7 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ...
    26 27 28 29 32 33 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53];

% get mechanical params for embryos and oocytes
% this requires calculating embryoUnkP and oocyteUnkP with
% "plotOocyteToEmbryoChange.m". We are missing 5,31,34
embryoUnkP = [embryoUnkP(1:4,:); NaN*ones(1,4); embryoUnkP(5:29,:); ...
    NaN*ones(1,4); embryoUnkP(30:31, :); NaN*ones(1,4); embryoUnkP(32:end,:)];
oocyteUnkP = [oocyteUnkP(1:4,:); NaN*ones(1,4); oocyteUnkP(5:29,:); ...
    NaN*ones(1,4); oocyteUnkP(30:31, :); NaN*ones(1,4); oocyteUnkP(32:end,:)];
embryoInfo.embryoUnkP = embryoUnkP;
embryoInfo.oocyteUnkP = oocyteUnkP;

% embryoUnkP = embryoInfo.embryoUnkP;
% oocyteUnkP = embryoInfo.oocyteUnkP;

for i = 1:length(embryoNums)
    % put oocyte and embryo mech params in E'i'
    embryoInfo.(['E' num2str(embryoNums(i))]).E1.oocyteParams = ...
        oocyteUnkP(mechNums(i),:);
    embryoInfo.(['E' num2str(embryoNums(i))]).E1.embryoParams = ...
        embryoUnkP(mechNums(i),:);    
end


% save(fileToSave, 'embryoInfo');

%% Plot amount of zona hardening against CG brightness after IVF


CGbright = [];
embryoK = [];
oocyteK = [];
zonaH = [];
colorMat = [];

for i = embryoNums
    
    currOP = embryoInfo.(['E' num2str(i)]).E1.oocyteParams;
    currEP = embryoInfo.(['E' num2str(i)]).E1.embryoParams; 
    
    CGbright = [CGbright embryoInfo.(['E' num2str(i)]).E1.brightParam];
    embryoK = [embryoK currEP(1)];
    oocyteK = [oocyteK currOP(1)];
    zonaH = [zonaH currEP(1)/currOP(1)];
    
end

figure(1);
clf;
set(gca, 'fontsize', 14);

% scatter3(zonaH, oocyteK, embryoK, 100, CGbright, 'filled');
scatter3(embryoK, oocyteK, CGbright, 100, [0 0 1], 'filled');
colormap default;
a = embryoNums';
b = num2str(a);
c = cellstr(b);
hold on;
text(embryoK + .1, oocyteK + .001, CGbright + .003, c);
xlabel('embryoK');
ylabel('oocyteK');
zlabel('CGbright');
view(0,90);


%% visualize CG images


figure(2);
clf;

for i = 1:length(embryoNums)
    
    subplot(7,7,i);
    hold on;
    imshow(embryoInfo.(['E' num2str(embryoNums(i))]).E1.ROI);
    title(['E' num2str(embryoNums(i))]);
%     colorbar;
    caxis([0 .01]);
    
end





















