% system identification applied to mouse embryo aspiration curves
% Applies 5 different models to the data:
% Maxwell (2 elements), Kelvin-Voigt (2 elements), SLS (3 elements),
% modified SLS (our model, 4 elements), and Wiechert (5 elements)
%
% Fit each model to each embryo, and calculate its error compared to all
% other embryos
%
% First: run the script plotAllMouseEmbryosSoFar
% This will result in 2 matrices aspirationDataN and aspirationDataV
% containing the data for all embryos

fileToSave = 'functions\TestModels\fitDataViableOnly.mat';
aspirationDataAll = aspirationDataV; %[aspirationDataV; aspirationDataN];
s = size(aspirationDataAll,1);
F = .325 * 6894.75729 * pi * (20 * 10^-6)^2; % pressure * area
plotInput = 1;

fMaxwell = zeros(1,s);
fKelvinVoigt = zeros(1,s);
fSLS = zeros(1,s);
fMod = zeros(1,s);
fFive = zeros(1,s);

% validation error (for each embryo on rest of embryos)
vMaxwell = zeros(s,s);
vKelvinVoigt = zeros(s,s);
vSLS = zeros(s,s);
vMod = zeros(s,s);
vFive = zeros(s,s);


%% 1. Maxwell model


for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));
    
    startParams = zeros(1,2);
    startParams(1) = .15; % k
    startParams(2) = .5; % n
    
    [paramsOut, fval] = MaxwellOptimizeParams(tData, xData, startParams, F, plotInput);
    fMaxwell(i) = fval;
    
    k = paramsOut(1);
    n = paramsOut(2);
    
    params = [F, n];
    options = odeset();
    
    [tOut, xOut] = ode45( @(t,x)MaxwellTest(t,x,params), tData, ...
        F/k, options);
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));
    
    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end
    
    % calculate fval for these params on all other embryos
    for j = 1:s
        xData = aspirationDataAll(j,:);
        vMaxwell(i,j) = sum(weighting.*(xOut' - xData).^2);
    end
end

%% 2. Kelvin-Voigt model


for i = 1:s
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));
    
    startParams = zeros(1,2);
    startParams(1) = .15; % k
    startParams(2) = .002; % n
    
    [paramsOut, fval] = KelvinVoigtOptimizeParams(tData, xData, startParams, F, plotInput);
    fKelvinVoigt(i) = fval;
    
    k = paramsOut(1);
    n = paramsOut(2);
    
    params = [F, k, n];
    options = odeset();
    [tOut, xOut] = ode45( @(t,x)KelvinVoigtTest(t,x,params), tData, ...
        0, options);
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));
    
    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end
    
    % calculate fval for these params on all other embryos
    for j = 1:s
        xData = aspirationDataAll(j,:);
        vKelvinVoigt(i,j) = sum(weighting.*(xOut' - xData).^2);
    end
    
end

%% 3. SLS model


for i = 1:size(aspirationDataAll,1)
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));
    
    startParams = zeros(3,1);
    startParams(1) = .05; % k0
    startParams(2) = .15; % k1
    startParams(3) = .005; % n
    
    [paramsOut, fval] = SLSOptimizeParams(tData, xData, startParams, F, plotInput);
    fSLS(i) = fval;
    
    k0 = paramsOut(1);
    k1 = paramsOut(2);
    n = paramsOut(3);
    
    params = [F, k0, k1, n];
    [tOut, xOut] = ode45( @(t,x)SLSTest(t,x,params), tData, ...
        F/(k0 + k1), options);
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));
    
    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end
        
    % calculate fval for these params on all other embryos
    for j = 1:s
        xData = aspirationDataAll(j,:);
        vSLS(i,j) = sum(weighting.*(xOut' - xData).^2);
    end

end

%% 4. Modified SLS (our model)


for i = 1:size(aspirationDataAll,1)
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));
    
    startParams = zeros(4,1);
    startParams(1) = .02; % k0
    startParams(2) = .15; % k1'
    startParams(3) = .001; % n0
    startParams(4) = 10; % n1
    
    [paramsOut, fval] = ModifiedSLSOptimizeParams(tData, xData, startParams, F, plotInput);
    fMod(i) = fval;
    
    k0 = paramsOut(1);
    k1 = paramsOut(2);
    n0 = paramsOut(3);
    n1 = paramsOut(4);
    
    params = [F, k0, k1, n0, n1];
    options = odeset();
    x0 = F/(k0 + k1);
    
    [tOut, xOut] = ode45( @(t,x)ModifiedSLSTest(t,x,params), tData, ...
        [0 x0 0 x0], options);
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));
    
    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end

    % calculate fval for these params on all other embryos
    for j = 1:s
        xData = aspirationDataAll(j,:);
        vMod(i,j) = sum(weighting.*(xOut(:,4)' - xData).^2);
    end

end

%% 5. Five-param model (Wiechert)

for i = 1:size(aspirationDataAll,1)
    
    i 
    xData = aspirationDataAll(i,:);
    tData = 0:.013:(.013*(length(xData) - 1));
    
    startParams = zeros(5,1);
    startParams(1) = .15; % k0
    startParams(2) = .03; % k1
    startParams(3) = .001; % k2
    startParams(4) = .5; % n0
    startParams(5) = .001; % n1
    
    [paramsOut, fval] = FiveParamOptimizeParams(tData, xData, startParams, F, plotInput);
    fFive(i) = fval;
    
    k0 = paramsOut(1);
    k1 = paramsOut(2);
    k2 = paramsOut(3);
    n0 = paramsOut(4);
    n1 = paramsOut(5);
    
    params = [F, k0, k1, k2, n0, n1];
    options = odeset();
    x0 = F/(k0 + k1 + k2);
    
    [tOut, xOut] = ode45( @(t,x)FiveParamTest(t,x,params), tData, ...
        [x0 0 x0 0 x0], options);
    
    % place extra weight on good fit in first data points
    weighting = 1*ones(1, length(xData));
    
    if length(xData) > 20
        weighting(1:10) = 2;
        weighting(1) = 10;
    else
        weighting(1:length(xData)-2) = 1;
    end
        
    % calculate fval for these params on all other embryos
    for j = 1:s
        xData = aspirationDataAll(j,:);
        vFive(i,j) = sum(weighting.*(xOut(:,5)' - xData).^2);
    end
    
end

%% Evaluate fitting error

mean(fMaxwell)
mean(fKelvinVoigt)
mean(fSLS)
mean(fMod)
mean(fFive)


figure(2);
clf;
h1 = scatter([2 2 3 4 5], log([mean(fMaxwell), mean(fKelvinVoigt), ...
    mean(fSLS), mean(fMod), mean(fFive)]), 150, [0 0 1], 'filled');
set(h1, 'Marker', 'd');
set(gca, 'XTick', [2 3 4 5]);
set(gca, 'XTickLabel', {'2', '3', '4', '5'});
set(gca, 'fontsize', 14);
title('Comparison of Different Bulk Models');
xlabel('number of model parameters');
ylabel('log(error)');
xlim([1.9 5.1]);
grid on;

% Evaluate validation error
vMeans = zeros(1,5);
vMeans(1) = mean(mean(vMaxwell));
vMeans(2) = mean(mean(vKelvinVoigt));
vMeans(3) = mean(mean(vSLS));
vMeans(4) = mean(mean(vMod));
vMeans(5) = mean(mean(vFive))

figure(2);
hold on;
h2 = scatter([2 2 3 4 5], log(vMeans), 150, [1 0 0]);
set(h2, 'Marker', 'o');
set(gca, 'XTick', [2 3 4 5]);
set(gca, 'XTickLabel', {'2', '3', '4', '5'});
set(gca, 'fontsize', 14);
title('Comparison of Different Bulk Models');
xlabel('number of model parameters');
ylabel('log(error)');
xlim([1.9 5.1]);
ylim([-28 -15])
grid on;

save(fileToSave, 'fMaxwell', 'fKelvinVoigt', 'fSLS', 'fMod', 'fFive', 'vMaxwell', 'vKelvinVoigt', 'vSLS', 'vMod', 'vFive');
legend([h1 h2], {'Fit error', 'Validation error'})














