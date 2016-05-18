% load mouse oocyte data

function outputStruct = loadDataOocytes()

outputStruct = struct;
numExperiments = 9;
experimentType = zeros(1,numExperiments); % 0 is non-fert, 1 is fert
dateList = cell(1,numExperiments);
numOocytes = zeros(1,numExperiments); % total num in experimental group
oocyteNums = cell(1,numExperiments); % actual labels of meas oocytes
fertInfo = cell(1,numExperiments);
blastForm = cell(1,numExperiments);
hatchInfo = cell(1,numExperiments);
maturationEnv = cell(1,numExperiments); % 0 = in vivo, 1 = KSOM, 2 = MM
morphologyInfo = cell(1,numExperiments); % 0 = GV, 1 = M1, 2 = MII, -1 = fragmented

% how much to add to k1 relative to what it would be if 0.2 psi were used
% for example, for .1 psi meas pressure, add .0314 to k1
k1ScaleFactor = zeros(1,numExperiments); 


% 3-21-16
% oocyte collection at 12 hrs post hCG
% mech measurement and fertilization at 16 hrs
% 0.1 psi pressure used
% 10 oocytes total in measured group (E1-10 were also measured post
% fertilization which I think harmed them; E11-20 were only measured before
% fertilization)
% controls were group cultured
% controls: 33/42 fertilized, 10/42 blast formation
experimentType(1) = 1;
dateList{1} = '3-21-16';
numOocytes(1) = 10; 
oocyteNums{1} = [11:20];
fertInfo{1} = [0 1 1 1 0 1 1 1 1 0];
blastForm{1} = [0 1 0 1 0 0 0 1 1 0];
hatchInfo{1} = [0 0 0 1 0 0 0 0 1 0];
maturationEnv{1} = zeros(1,10);
k1ScaleFactor(1) = .0314;
morphologyInfo{1} = [1 2 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];


% 3-28-16
% oocyte collection at 12 hrs post hCG
% mech measurement and fertilization at 16 hrs
% 0.1 psi pressure used
% 24 oocytes total in measured group, only measured before fertilization
% controls were group cultured
% controls: 24/40 fertilized, 20/40 blast formation
experimentType(2) = 1;
dateList{2} = '3-28-16';
numOocytes(2) = 24; 
oocyteNums{2} = 1:24;
fertInfo{2} = [zeros(1,5) 1 0 0 0 1 0 0 1 1 0 0 1 1 0 0 0 1 1 1];
blastForm{2} = [zeros(1,12) 1 1 0 0 0 1 0 0 0 1 1 1];
hatchInfo{2} = [zeros(1,12) 1 1 0 0 0 1 0 0 0 1 1 1];
maturationEnv{2} = zeros(1,24);
k1ScaleFactor(2) = .0314;
morphologyInfo{2} = [2 2 2 2 2 2 2 1 2 2 1 2 2 1 2 2 2 2 1 2 2 2 1 2];


% 10-8-15 
% In vitro maturation only
% 0.2 psi used
% mech group was measured at 0,8,16,24 hrs (so only take _24 expt), where 0
% is 5 hrs post hCG.
% morphology is measured at 24 hr mark
% exp group was measured only once at 24 hrs, but brought outside incubator
%   while mech group was being measured
% inc group was kept inside incubator and only measured once at 24 hrs
experimentType(3) = 0;
dateList{3} = '10-8-15';
numOocytes(3) = 65; 
oocyteNums{3} = 1:65;
fertInfo{3} = NaN*ones(1,65);
blastForm{3} = NaN*ones(1,65);
hatchInfo{3} = NaN*ones(1,65);
maturationEnv{3} = [1 1 1 1 1 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 ...
    ones(1,16) 2*ones(1,12) 2*ones(1,9) ones(1,8)];
k1ScaleFactor(3) = 0;
morphologyInfo{3} = [2 0 -1 2 2 2 2 2 2 1 1 -1 2 2 1 1 2 1 1 2 ...
    1 1 0 0 1 1 0 1 1 2 0 1 1 2 2 2 2 2 2 1 2 -1 2 2 1 2 2 2 ...
    2 2 2 1 1 1 0 1 1 2 2 2 1 2 1 1 -1];


% 10-21-15 
% In vitro maturation only
% 0.2 psi used
% morphology is measured at 24 hr mark
% mech group was measured at 0,8,16,24 hrs (so only take _24 expt), where 0
% is 5 hrs post hCG.
% exp group was measured only once at 24 hrs, but brought outside incubator
%   while mech group was being measured
% inc group was kept inside incubator and only measured once at 24 hrs
experimentType(4) = 0;
dateList{4} = '10-21-15';
numOocytes(4) = 71; 
oocyteNums{4} = 1:71;
fertInfo{4} = NaN*ones(1,71);
blastForm{4} = NaN*ones(1,71);
hatchInfo{4} = NaN*ones(1,71);
maturationEnv{4} = [ones(1,10) 2*ones(1,20) ones(1,11) ...
    2*ones(1,10) 2*ones(1,10) ones(1,10)];
k1ScaleFactor(4) = 0;
morphologyInfo{4} = [0 2 -1 1 1 0 2 1 -1 2 2 0 2 2 1 2 0 2 2 2 ...
    2 0 -1 2 2 2 2 2 2 2 2 2 1 1 2 2 1 1 2 2 2 1 0 2 2 2 2 1 2 2 2 1 ...
    2 2 1 2 2 1 -1 2 1 1 1 1 1 2 -1 0 2 1 1];


% 12-2-15
% attempted IVF but failed
% 0.2 psi used
% oocytes matured in vivo and collected 16 hrs post hCG
% mech measurement done 16 hrs post hCG
experimentType(5) = 0;
dateList{5} = '12-2-15';
numOocytes(5) = 30; 
oocyteNums{5} = 1:30;
fertInfo{5} = zeros(1,30);
blastForm{5} = zeros(1,30);
hatchInfo{5} = zeros(1,30);
maturationEnv{5} = zeros(1,30);
k1ScaleFactor(5) = 0;
morphologyInfo{5} = [1 2 0 1 1 2 2 2 2 1 1 2 1 2 2 2 1 1 2 1 2 2 1 1 ...
    1 1 2 2 2 1];


% 12-10-15
% attempted IVF but failed
% 0.2 psi used
% morphology is measured at 16 hr mark
% oocytes collected 3 hrs post hCG, matured in MM (1-20) or KSOM (21-40)
% mech measurement and fertilization done at 16 hrs post hCG
experimentType(6) = 0;
dateList{6} = '12-10-15';
numOocytes(6) = 40; 
oocyteNums{6} = 1:40;
fertInfo{6} = NaN*ones(1,40);
blastForm{6} = NaN*ones(1,40);
hatchInfo{6} = NaN*ones(1,40);
maturationEnv{6} = [2*ones(1,20) ones(1,20)];
k1ScaleFactor(6) = 0;
morphologyInfo{6} = [2 2 1 0 1 2 1 1 2 1 2 1 1 2 2 2 2 1 2 0 1 1 0 1 ...
    1 0 1 1 1 2 1 0 1 0 1 1 1 1 1 1];


% 1-13-16
% attempted IVF but failed
% 0.2 psi used
% oocytes collected 3 hrs post hCG, matured in MM (1-10) or KSOM (11-20)
% filenames ending in _16 correspond to mech measurement at 16 hrs post hCG
% mech measurement done at 5,16 hrs post hCG, fertilization done at 16 hrs post hCG
experimentType(7) = 0;
dateList{7} = '1-13-16';
numOocytes(7) = 20; 
oocyteNums{7} = 1:20;
fertInfo{7} = NaN*ones(1,20);
blastForm{7} = NaN*ones(1,20);
hatchInfo{7} = NaN*ones(1,20);
maturationEnv{7} = [2*ones(1,10) ones(1,10)]; % this is just a guess - these may be switched
k1ScaleFactor(7) = 0;
morphologyInfo{7} = [0 1 1 2 1 2 1 0 2 1 1 2 0 0 2 1 2 2 1 2];


% 1-20-16
% attempted IVF but failed
% 0.2 psi used
% oocytes collected 16 hrs post hCG
% mech measurement done at 16 hrs post hCG, fertilization done at 16 hrs post hCG
experimentType(8) = 0;
dateList{8} = '1-20-16';
numOocytes(8) = 20; 
oocyteNums{8} = 1:20;
fertInfo{8} = zeros(1,20);
blastForm{8} = zeros(1,20);
hatchInfo{8} = zeros(1,20);
maturationEnv{8} = zeros(1,20);
k1ScaleFactor(8) = 0;
morphologyInfo{8} = [1 2 1 2 2 2 1 2 2 1 1 2 1 2 2 1 2 1 1 2];


% 2-2-16
% attempted IVF, only worked on controls
% 0.2 psi used
% oocytes collected 12 hrs post hCG
% mech measurement done at 14 hrs post hCG, fertilization done at 16 hrs post hCG
experimentType(9) = 0;
dateList{9} = '2-2-16';
numOocytes(9) = 15; 
oocyteNums{9} = 1:15;
fertInfo{9} = NaN*ones(1,15);
blastForm{9} = NaN*ones(1,15);
hatchInfo{9} = NaN*ones(1,15);
maturationEnv{9} = zeros(1,15);
k1ScaleFactor(9) = 0;
morphologyInfo{9} = [1 1 2 2 2 2 2 2 2 2 1 1 2 2 2];



outputStruct.experimentType = experimentType;
outputStruct.numExperiments = numExperiments;
outputStruct.dateList = dateList;
outputStruct.numOocytes = numOocytes;
outputStruct.oocyteNums = oocyteNums;
outputStruct.fertInfo = fertInfo;
outputStruct.blastForm = blastForm;
outputStruct.hatchInfo = hatchInfo;
outputStruct.maturationEnv = maturationEnv;
outputStruct.k1ScaleFactor = k1ScaleFactor;
outputStruct.morphologyInfo = morphologyInfo;


