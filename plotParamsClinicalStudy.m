% plot params for all patients in clinical study
% Livia Zarnescu Yanez 11-10-15

params = loadDataClinicalStudy();
numParticipants = params.numParticipants;
participantIDs = params.participantIDs;
numEmbryos = params.numEmbryos;
outcomeInfo = params.outcomeInfo;
d3M = params.d3M;
blastM = params.blastM;
ICSI = params.ICSI;
PGD = params.PGD;
gender = params.gender;
totalNumEmbryos = sum(numEmbryos);
procDataPath = 'C:\Users\Livia\Desktop\IVF\Processed Data\Human\';

% init param vectors
mList = [];
k0list = [];
k1list = [];
n0list = [];
n1list = [];
taulist = [];
gList = [];
colorMat = [];
currE = 0;

ptsToPlot = [1 1 1 1 1 1 1 1];
aText = [];
cText = {};

% load in params by participant and embryo num
for i = 1:numParticipants
    if ptsToPlot(i)
        for j = 1:numEmbryos(i)
            
            currE = currE + 1;
            aText = [aText currE];
            cText = {cText{:}, ['P', num2str(i), '\_E', num2str(j)]};
            
            currDataPath = [procDataPath participantIDs{i} '\' ...
                participantIDs{i} '_E' num2str(j) '.mat'];
            
            % save embryo params and color
            if exist(currDataPath, 'file')
                load(currDataPath);
                mList = [mList outcomeInfo{i}(j)];
                k0list = [k0list k0];
                k1list = [k1list k1];
                n0list = [n0list n0];
                n1list = [n1list n1];
                taulist = [taulist tau];
                gList = [gList gender{i}(j)];
                
                if outcomeInfo{i}(j)
                    
                    colorMat = [colorMat; [0 .6 0]];
%                     if gender{i}(j) == 1
%                         colorMat = [colorMat; [0 .6 0]];
%                     elseif gender{i}(j) == 0
%                         colorMat = [colorMat; [.9 .4 0]];
%                     else
%                         colorMat = [colorMat; [.6 .6 0]];
%                     end
                    %                     if ICSI{i}(j) == 1
                    %                         colorMat = [colorMat; [0 .6 0]];
                    %                     else
                    %                         colorMat = [colorMat; [.6 .6 0]];
                    %                     end
                else
                    %                     if ICSI{i}(j) == 1
                    %                         colorMat = [colorMat; [0 0 .6]];
                    %                     else
                    %                         colorMat = [colorMat; [0 .6 .6]];
                    %                     end
                    colorMat = [colorMat; [0 0 .6]];
                end
                
            else
                mList = [mList NaN];
                k0list = [k0list NaN];
                k1list = [k1list NaN];
                n0list = [n0list NaN];
                n1list = [n1list NaN];
                taulist = [taulist NaN];
            end
        end
    end
end


% now plot!
p1 = k1list;
p2 = n1list;
p3 = taulist;

figure(1); clf;
% h = scatter3(p1, p2, p3, 200, colorMat, 'filled');
h = scatter(p1, p2, 200, colorMat, 'filled');

set(h, 'Marker', 'o');
set(gca, 'FontSize', 14);
title('Mechanical parameters predict blastocyst formation');
xlabel('k1 parameter');
ylabel('n1 parameter');
% zlabel('tau parameter');
% axis([min(p1) max(p1) min(p2) max(p2) min(p3) max(p3)]);
% set(gca, 'yscale', 'log')
axis([min(p1) max(p1) min(p2) max(p2)]);
set(gca, 'zscale', 'linear');
grid on;

% b = num2str(aText');
% c = cellstr(b);
% dx = -0.001; dy = 0.015;% dz = .003; % displacement so the text does not overlay the data points
% hold on;
% text(p1+dx, p2+dy, cText'); %, p3+dz, c);







