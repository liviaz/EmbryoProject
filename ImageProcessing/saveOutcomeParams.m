% Save outcome params to .mat file
%
%

function [] = saveOutcomeParams()

clear all;

dishNum = [ones(1,11) 2*ones(1,8)];
wellNum = {'B1', 'D1', 'E1', 'A2', 'C2', 'E2', 'B3', 'D3', 'A4', 'C4', 'E4', ...
    'B1', 'D1', 'E1', 'A2', 'C2', 'E2', 'B3', 'D3'};

blastForm = [1 0 0 1 1 0 1 1 0 0 1 ...
             0 0 0 0 0 1 0 0];
blastMorph = [9 0 0 3 10 0 3 8 0 0 3 ...
             0 0 0 0 0 7 0 0];

chromNorm = [1 0 0 1 1 0 NaN 1 1 1 1 NaN 0 NaN 1 0 1 NaN 0];


save('C:\Users\Livia\Desktop\Human embryo map\OutcomeParams.mat');