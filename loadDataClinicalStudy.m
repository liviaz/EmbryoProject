% load outcome data clinical study

function outputStruct = loadDataClinicalStudy()

outputStruct = struct;
numParticipants = 8;
participantIDs = cell(1,numParticipants);
numEmbryos = zeros(1,numParticipants);
outcomeInfo = cell(1,numParticipants);
d3M = cell(1,numParticipants);
blastM = cell(1,numParticipants);
ICSI = cell(1,numParticipants);
PGD = cell(1,numParticipants);
gender = cell(1,numParticipants); % 0 = male, 1 = female

% EMECH-001
% 10-3-15
% E1 and E2 were non-ferts
participantIDs{1} = 'MECH001';
numEmbryos(1) = 9; 
outcomeInfo{1} = [0 0 1 1 1 0 1 0 0];
d3M{1} = {}; % don't have for this patient
blastM{1} = {'', '', '6BB', '5AA', '5BB', '', '5BB', '', ''};
ICSI{1} = ones(1,9);
PGD{1} = NaN*ones(1,9);
gender{1} = NaN*ones(1,9);

% EMECH-002
% exited from study due to low oocyte count
participantIDs{2} = 'MECH002';
numEmbryos(2) = 0; 
outcomeInfo{2} = [];
d3M{2} = {};
blastM{2} = {};
ICSI{2} = [];
PGD{2} = [];
gender{2} = [];

% EMECH-003
% 10-23-15
% all blasts frozen for PGD and future transfer
participantIDs{3} = 'MECH003';
numEmbryos(3) = 14; 
outcomeInfo{3} = [1 1 1 1 1 1 1 1 1 0 0 1 0 0];
d3M{3} = {'8I', '7I', '8I', '8II', '8II', '8I', '8I', '5III', '5II', ...
    '5III', '4III', '8II', '7III', '8I'};
blastM{3} = {'5BB', '5BA', '3BB', '3BB', '4BA', '5AA', '5BB', '5BB', ...
    '5BB', '', '2CC', '5AB', '5CB', ''};
ICSI{3} = ones(1,14);
PGD{3} = [1 1 1 1 1 1 0 0 1 NaN NaN 1 NaN NaN];
gender{3} = [0 1 0 0 1 1 1 1 0 NaN NaN 0 NaN NaN];

% EMECH-004
% 10-23-15
% all discarded due to poor quality
participantIDs{4} = 'MECH004';
numEmbryos(4) = 13;
outcomeInfo{4} = [zeros(1,11) 1 0];
d3M{4} = {'8III', '6III', '8II', '8II', '8I', '4IV', '5III', '6III', ...
    '4IV', '8II', '5III', '8II', '7III'};
blastM{4} = {'2CC', '', '1CC', '1', '', '', '', '', '', '', '', '5CC', ''};
ICSI{4} = ones(1,13);
PGD{4} = NaN*ones(1,13);
gender{4} = NaN*ones(1,13);

% EMECH-005
% 10-25-15
participantIDs{5} = 'MECH005';
numEmbryos(5) = 10;
outcomeInfo{5} = [1 0 1 1 0 1 0 1 0 1];
d3M{5} = {'8II', '4I', '7I', '7I', '5V', '7II', '7II', '8II', '4V', '7II'};
blastM{5} = {'5BB', '', '6BB', '6AA', '', '5BB', '', '5BB', '', '5BB'};
ICSI{5} = [0 0 0 ones(1,7)];
PGD{5} = NaN*ones(1,10);
gender{5} = NaN*ones(1,10);

% EMECH-006
% 11-9-15, last 2 embryos were 1PNs
participantIDs{6} = 'MECH006';
numEmbryos(6) = 7;
outcomeInfo{6} = [0 0 1 1 0 1 0];
d3M{6} = {'4I', '8III', '6II', '8II', '9III', '7II', '4IV'};
blastM{6} = {'', '', '6BB', '5BC', '', '3BC', ''};
ICSI{6} = [ones(1,7)];
PGD{6} = [NaN NaN 1 1 NaN 1 NaN];
gender{6} = [NaN NaN 1 1 NaN 1 NaN];

% EMECH-007
% has not had retrieval yet
participantIDs{7} = 'MECH007';
numEmbryos(7) = 0;
outcomeInfo{7} = [];
d3M{7} = {};
blastM{7} = {};
ICSI{7} = [];
PGD{7} = [];
gender{7} = [];

% EMECH-008
% 11-11-15
participantIDs{8} = 'MECH008';
numEmbryos(8) = 11;
outcomeInfo{8} = [0 0 0 0 1 0 0 0 1 0 0];
d3M{8} = {'3II', '3II', '8IV', '8II', '9II', '7II', '4V', '8II', ...
    '8II', '7IV', '8I'};
blastM{8} = {'', '', '', '', '6BB', '', '', '', '1', '', ''};
ICSI{8} = zeros(1,11);
PGD{8} = NaN*ones(1,11);
gender{8} = NaN*ones(1,11);

outputStruct.numParticipants = numParticipants;
outputStruct.participantIDs = participantIDs;
outputStruct.numEmbryos = numEmbryos;
outputStruct.outcomeInfo = outcomeInfo;
outputStruct.d3M = d3M;
outputStruct.blastM = blastM;
outputStruct.ICSI = ICSI;
outputStruct.PGD = PGD;
outputStruct.gender = gender;






