% Process and Analyze Multi-Zone Visco-NPS Signals
% Andre Lai | June 20, 2022
%
% Edit parameters/sections accordingly. 
% Execute code section by section, in order.
% Execute each section only once.

%% Global Fixed Parameters For Rheological Analysis

clc; clearvars;

% *** EDIT the following parameters: ***
L = 26550; % full pore-to-pore channel length [um]
Lzone = [6, 6, 6, 6]; % length of each contraction zone [mm]
Lpore = 1; % length of the sizing pore [mm]
numZones = 5; % pore + number of sinusoidal contraction zones (6 is max)
periods = [4, 8, 12, 16] * 3; % number of periods in each zone [unitless]

% ----------- create dictionary of dP values ------------------------------
% dP is the COMSOL modeled average pressure diff. around cell [Pa]
% dP is dependent on channel geometry, input pressure, and cell size
keySet = {'W9.25-1PSI','W9.25-2PSI','W9.25-3PSI','W9.25-4PSI',...
          'W11.25-1PSI','W11.25-2PSI','W11.25-3PSI','W11.25-4PSI'};
f9_1 = @(x) 2.32*exp(0.15*x) + (4.6e-5)*exp(0.67*x);
f9_2 = @(x) 4.64*exp(0.15*x) + (9.2e-5)*exp(0.67*x);
f9_3 = @(x) 6.62*exp(0.15*x) + (1.3e-4)*exp(0.67*x);
f9_4 = @(x) 8.94*exp(0.15*x) + (1.8e-4)*exp(0.67*x);
f11_1 = @(x)  5.09*exp(0.08*x) + (5.5e-3)*exp(0.39*x);
f11_2 = @(x) 10.16*exp(0.08*x) + (1.1e-2)*exp(0.39*x);
f11_3 = @(x) 14.50*exp(0.08*x) + (1.5e-2)*exp(0.39*x);
f11_4 = @(x) 19.54*exp(0.08*x) + (2.1e-2)*exp(0.39*x);
valueSet = {f9_1, f9_2, f9_3, f9_4, f11_1, f11_2, f11_3, f11_4};
dPDictionary = containers.Map(keySet,valueSet);

% ---------- create dictionary of Deff values -----------------------------
% Deff is the % experimentally derived effective diameter of the pore [um]
% Values are calibrated to the avg. cell diameters as measured by the
%   Millapore Scepter (2.0)
% Values can be obtained by running the EffectiveDiameter.m script
keySet = {'W9.25', 'W11.25'};
valueSet = [16.65, 19.25];
DeffDictionary = containers.Map(keySet,valueSet);

% ---------- create dictionary of device geometries -----------------------
% arrays in valueSet stored as [w_max, w_min, w_np]
% w_max: maximum channel width [um]
% w_min: minimum channel width [um]
% w_np: channel width at node-pore [um]
keySet = {'W9.25', 'W11.25'};
valueSet = {[11.5,7,17],[14,8.5,21]};
wDictionary = containers.Map(keySet,valueSet);

clear keySet valueSet f9_1 f9_2 f9_3 f9_4 f11_1 f11_2 f11_3 f11_4

%% Load .mat Raw Data File
% "data" matrix from [.mat] file should be organized as a [n x 1] vector in
% the order: [Current, pore, V1, V2, V3, and so on ...] for each sample.

% *** EDIT [filreDir] & [fileName] to select desired data file ***
fileDir = '';
fileName = '';

fileTestName = fileName(1:5);

% extract experimental date
date = fileDir(70:77); % *EDIT* if fileDir pathway format has changed

% extract device ID
deviceID = fileDir(81:88); % *EDIT* if fileDir pathway format has changed

% extract cell type
possibleCellTypes = ["MCF7", "MDA", "MCF10A"];

for type = possibleCellTypes
    if contains(fileName, type)
        cellType = type;
    end
end

cellType = fileName(7:13);

% extract operating pressure for this data file
possiblePressureValues = ["1PSI", "2PSI", "3PSI", "4PSI"];

for pressure = possiblePressureValues
    if contains(fileName, pressure)
        inputPressure = char(pressure);
    end
end

% extract channel width for this device 
possibleWidthValues = ["W9.25", "W11.25"];

for width = possibleWidthValues
    if contains(fileDir, width)
        deviceWidth = char(width);
    end
end

deviceWidth = 'W9.25';

% determine Deff from deviceWidth
Deff = DeffDictionary(deviceWidth);

% determine Deff from deviceWidth
w = wDictionary(deviceWidth);
w_max = w(1); % maximum channel width [um]
w_min = w(2); % minimum channel width [um]
w_np = w(3); % channel width at node-pore [um]

% determine dP from extracted inputPressure and deviceWidth
dP = dPDictionary(strcat(deviceWidth, "-", inputPressure));

% load the data 
load(strcat(fileDir, fileName))

fprintf('Data Loaded: %s - Device %s - %s\n',date, deviceID, fileName)

clear demodSignal freqVin Vpp possibleCellTypes type possiblePressureValues ...
      pressure possibleWidthValues width dPDictionary DeffDictionary w ...
      wDictionary DhydDictionary

%% Retrieve Data [V,I] from .mat Filter and Calculate Resistances [R, Rf]

% 1. Retrieve Raw Current Signal ------------------------------------------

% multiple raw voltage data with ampsPerVolt parameter to get current data
I = transpose(data(1:numZones+1:end) * ampsPerVolt);

% 2. Retrieve Raw Voltage Signals -----------------------------------------

% intialize raw voltage matrix
V = zeros(length(data)/(numZones+1), numZones);

% populate matrix with raw voltage signal for pore and each zone
for i = 1:numZones
    V(:,i) = transpose(data(i+1:numZones+1:end));
end

% create total voltage vector
Vtotal = sum(V(:,1:numZones),2);

% *optional* trim signal **************************************************
if (0)
    start = 1; % *** EDIT Start
    stop = 1; % *** Edit End
    I = I(start:end);
    V = V(start:end,:);
    fprintf('Raw signal has been trimmed. Double check the signal.\n')
end

% 3. Compute Time Vector At Original Sample Rate --------------------------

t = (0:size(I,1)-1).' ./ sampleRate; 
 
% 4. Calculate Resistance Signals -----------------------------------------

fpass = 200; % passband frequency of the lowpass filter (in Hz)
window = 15; % averaging window width (in units of samples) 
R = V ./ I;

% 5. Filter Resistance Signals --------------------------------------------

Rf = smoothdata(lowpass(R, fpass, sampleRate), 'movmean', window);
%Rf_band = smoothdata(bandpass(R, [10 fpass], sampleRate), 'movmean', window);

clear ampsPerVolt data i fpass window

%% Identify Baseline Resistances
% *note* running this section takes approximately ~ 2 min

tic % start timer

% baseline parameters
bias = 0.0;
envelopeWindow = round(sampleRate/50);
window = sampleRate*3;

% initialize baseline resistance matrix
Rbase = zeros(length(R), numZones);

% find baselines
for i = 1:numZones
    
    % display progress in command window
    if i == 1
        lineLength = fprintf('Performing Baseline Calculation On Zone %d of %d...\n', i, numZones);
    else
        fprintf(repmat('\b',1,lineLength))
        lineLength = fprintf('Performing Baseline Calculation On Zone %d of %d...\n', i, numZones);
    end
        
    % find the upper envelope of the signal
    [upper,lower] = envelope(Rf(:,i), envelopeWindow, 'peak');
    
    % calculate baseline from signal envelopes and shift the baseline by the
    % [bias] percent of the range between the upper envelope median value and
    % the baseline median value
    baseline = movmedian((upper+lower)/2, window);
    baseline = baseline + bias * (median(upper) - median(baseline));
    Rbase(:,i) = baseline;
end

% subtract the baseline
Rnorm = Rf-Rbase;

fprintf(repmat('\b',1,3))
fprintf("\nBaseline Identification Complete - ")
toc % end timer

clear i bias envelopeWindow window upper lower baseline lineLength

%% *Optional* Plot Baseline Resistance

% plot baselines ----------------------------------------------------------
dF = 4; % downsampling factor, downsampling improves plotting time & performance
figure(1)
tplot = tiledlayout(numZones,1);
tplot.TileSpacing = 'compact';
tplot.Padding = 'compact';
set(gcf,'Position',[50 400 500 800])
colors = [[0.3, 0.3, 0.3]; [0, 0.4470, 0.7410];...
          [0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];....
          [0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880]];

for i = 1:numZones
    
    nexttile;
    plot(downsample(t,dF), downsample(R(:,i),dF), 'Color', [160/255,160/255,160/255], 'linestyle', ':'); hold on
    plot(downsample(t,dF), downsample(Rf(:,i),dF), 'Color', colors(i,:), 'linewidth', 2); hold on
    plot(downsample(t,dF), downsample(Rbase(:,i),dF), 'Color', 'red', 'linewidth', 1); hold on
    
end

titleTop = sprintf("Resistance Signals & Baseline\n");
titleBottom = sprintf("Cell Type: %s | Device: %s %s | %s", cellType, deviceID, deviceWidth, inputPressure);
titleText = strcat(titleTop, titleBottom);
title(tplot, titleText, 'FontSize',20); 
xlabel(tplot,'Time [s]','FontSize',20)
ylabel(tplot,'Resistance [\Omega]', 'FontSize',20)

% link x-axes
allaxes = findobj(gcf, 'type','axes');
linkaxes(allaxes, 'x');

% misc. adjustments
axis tight;

% plot overlapping normalized resistances ---------------------------------
figure(2)
set(gcf,'Position',[550 700 500 300])

for i = 1:numZones
    plot(downsample(t,dF), downsample(Rnorm(:,i),dF)./1000,'Color',colors(i,:),'LineWidth',2)
    hold on
end

titleTop = sprintf("Normalized Resistance\n");
titleBottom = sprintf("Cell Type: %s | Device: %s %s | %s", cellType, deviceID, deviceWidth, inputPressure);
titleText = strcat(titleTop, titleBottom);
title(titleText, 'FontSize',20); 
xlabel('Time [s]','FontSize',20)
ylabel('Normalized Resistance [k\Omega]', 'FontSize',20)
legend('Pore', 'Zone 1', 'Zone 2', 'Zone 3', 'Zone 4', 'Zone 5','FontSize',14,'Location','Southeast');

% misc. adjustments
axis tight;
ylim([-100 200]);

% plot tiled normalized resistances ---------------------------------------
figure(3)
ttplot = tiledlayout(numZones,1);
ttplot.TileSpacing = 'compact';
ttplot.Padding = 'compact';
set(gcf,'Position',[1050 400 800 800])

for i = 1:numZones
    
    nexttile;
    plot(downsample(t,dF), downsample(Rnorm(:,i)./1000,dF),'Color',colors(i,:),'LineWidth',1.5)
    
    %xlim([450 490])
    %ylim([-4 4])
    if i ~= numZones
        xticks([])
    end
    %yticks([-4 -2 0 2 4])
    set(gca,'FontSize',18)
end

titleTop = sprintf("Normalized Resistance\n");
titleBottom = sprintf("Cell Type: %s | Device: %s %s | %s", cellType, deviceID, deviceWidth, inputPressure);
titleText = strcat(titleTop, titleBottom);
title(ttplot, titleText, 'FontSize',24); 
xlabel(ttplot,'Time [s]','FontSize',24)
ylabel(ttplot,'Normalized Resistance [k\Omega]', 'FontSize',24)

% link x-axes
allaxes = findobj(gcf, 'type','axes');
linkaxes(allaxes, 'x');
% misc. adjustments
axis tight;

clear titleTop titleBottom titleText colors tplot i allaxes ttplot

%% Identify Cell Transit Events via Baseline Intersect Intervals
% Using signal from Pore to identify cell transit events
% 
% Algorithm Approach: Determine where the signal intersects the calculated 
% baseline signal. Calculate the interval between the identified intersects. 
% Large intervals above the preset threshold are flagged as a cell transit 
% event. 

% *EDIT* Tunable parameters to adjust (False Positve / False Negative)
%        cell transit event identifications.
% 1.) Interval Threshold: The minimum interval size to be flagged as an event. 
%     Value is: integer, the number of samples. Dependent on inlet pressure.
      intervalThreshold = 350; 
% 2.) Proximal Threshold: The minimum interval between flagged events.
%     Value is: integer, the number of samples. Dependent on inlet pressure. 
      proximalThreshold = 800; 

% baseline signal
baseline = Rbase(:,1);

% pore filtered resistance signal
Rpore = Rf(:,1);

% find baseline intersections
baseIntersects = false(length(Rpore)-1,1);
for i = 1:length(Rpore)-1
    if (Rpore(i) <= baseline(i) && Rpore(i+1) >= baseline(i)) || ...
       (Rpore(i) >= baseline(i) && Rpore(i+1) <= baseline(i))
        baseIntersects(i) = 1;
    end
end
baseIntersects = find(baseIntersects==1);

% calculate the intervals between the identified intersections
intervals = zeros(length(baseIntersects),1);
for i = 2:length(baseIntersects)
    intervals(i) = baseIntersects(i) - baseIntersects(i-1);
end

% find indicies and flag record those indices as events where the identified 
% intersection interval is greater than the [intevalThreshold]
events = zeros(length(baseIntersects),1);
for z = 2:length(intervals)
    if intervals(z) > intervalThreshold
        events(z-1) = baseIntersects(z-1);
    end
end
events = nonzeros(events);

% filter out events that are within the [proximalThreshold]
remove = false(length(events),1);
for i = 2:length(events)-1
    if events(i+1)-events(i) < proximalThreshold || ...
       events(i)-events(i-1) < proximalThreshold
        remove(i) = 1;
    end
    
end
events(remove) = [];

fprintf("%d total events identified.\n", length(events));

clear i j k z bias intervalThreshold proximalThreshold window upper lower ... 
      baseIntersects intervals remove Rpore baseline
  
%% *Optional* Plot Event Identification
% Plot identified cell transit events to review false negatives. Adjust
% tunable parameters in previous section accordingly. 

% initiate figure
figure(1)

% add plot data
plot(t, Rf(:,1)); 
hold on
plot(t, Rbase(:,1), 'linewidth', 2, 'linestyle', ':'); 
hold on
scatter(t(events), Rf(events,1), 44, 'filled');

% add title
titleTop = sprintf("Event Identifications\n");
titleBottom = sprintf("Cell Type: %s | Device: %s %s | %s", ...
                       cellType, deviceID, deviceWidth, inputPressure);
titleText = strcat(titleTop, titleBottom);
title(titleText, 'FontSize',16); 

% add axis labels
xlabel('Time [s]','FontSize',14); ylabel('Resistance [Ohms]','FontSize',14);

% add legend
legend('Pore','Baseline','Event','FontSize',14,'Location','Northeast');

% misc. adjustments
axis tight; grid on;

clear titleTop titleBottom titleText

%% Inititate GUI for Manual Filter of Auto Event Identification
% Continuously plot identified cell transit events to manually reject false positives. 

% number of samples to include in event signal before and after identified 
% cell transit event
keySet = {'1PSI', '2PSI', '3PSI', '4PSI'};
plotBufferValueSet = [1000, 600, 500, 250]*5;
plotBufferDict = containers.Map(keySet,plotBufferValueSet);

plotBuffer = plotBufferDict(inputPressure); 

% number of samples after the identified start time anticipated to be the 
% end of the cell transit event
keySet = {'1PSI', '2PSI', '3PSI', '4PSI'};
endBufferValueSet = [18000, 12000, 10000, 5000]*10;
endBufferDict = containers.Map(keySet,endBufferValueSet);

endBuffer = endBufferDict(inputPressure);

% initiate manual selection GUI
EventSelectionGUI(t, V(:,2), events, fileDir, fileTestName, plotBuffer, endBuffer);
              
clear keySet plotBufferValueSet plotBufferDict endBufferValueSet endBufferDict

%% Parse Out Resistance Signal for Each Selected Event

% retrieve selected events 
% loads a [nx1] matrix named [selectedEvents] into the workspace
load(strcat(fileDir, fileTestName, '-EventSelection.mat'))
fprintf("%d total events selected.\n", length(selectedEvents));

% initialize 3D matrix with format:
%{
 each column: deltaRpore deltaR1 deltaR2 deltaR3 and so on...
              baseRpore  baseR1  baseR2  baseR3  and so on... time
 each row: data sample
 each page: different selected event
%}
EventSignals = zeros(2*plotBuffer+endBuffer+1, numZones*2+1, length(selectedEvents));

% populate 3D matrix with detrened event signals
for eventNum = 1:length(selectedEvents)
    
    startSample = selectedEvents(eventNum) - plotBuffer;
    if startSample < 1
        startSample = 1;
    end
    endSample = startSample + 2*plotBuffer+endBuffer; 
    if endSample > length(R)
        endSample = length(R);
    end
    
    % first [numZones] columns in the [EventSignals] matrix are the delta resistance
    for z = 1:numZones
        % Not using Rnorm due to inconsistencies in a global baseline identification
        % EventSignals(:,z,eventNum) = Rnorm(startSample:endSample,z);
        
        % Using a local baseline adjustment instead
        EventSignals(:,z,eventNum) = detrend(Rf(startSample:endSample,z),1);
    end

    % *optional* additional baseline correction
    % further remove any trends by manually shifting begining and end median values to 0
    for z = 1:numZones
        offset = linspace(median(EventSignals(1:1000,z,eventNum)), ...
                      median(EventSignals(end-1000:end,z,eventNum)), ...
                      length(EventSignals(:,z,eventNum)));
        EventSignals(:,z,eventNum) = EventSignals(:,z,eventNum) - offset.';
    end
    
    % next [numZones] columns in the [EventSignals] matrix are the baseline 
    % resistance vectors for each zone
    for z = 1:numZones
        EventSignals(:,numZones+z,eventNum) = Rbase(startSample:endSample,z);
    end                              
    
    % last column in the [EventSignals] matrix is the time vector 
    EventSignals(:,numZones*2+1,eventNum) = t(startSample:endSample);
    
end

clear i z startSample endSample eventNum offset

%% *Optional* Plot Individual Events

figure(1)

tplot = tiledlayout(numZones,1);
tplot.TileSpacing = 'compact';
tplot.Padding = 'compact';

% plot parameters
eventNum = 3; % choose which event to plot
timeSignal = EventSignals(:,numZones*2+1,eventNum);
startTime = timeSignal(1); endTime = timeSignal(end);
scale = 10^3;
yStart = min(min(EventSignals(:,1:numZones,eventNum)))./scale;
yEnd = max(max(EventSignals(:,1:numZones,eventNum)))./scale;
colors = [[160/255,160/255,160/255];[0, 0.4470, 0.7410]; ...
          [0.8500, 0.3250, 0.0980] ;[0.9290, 0.6940, 0.1250]; ...
          [0.4940, 0.1840, 0.5560] ;[0.4660, 0.6740, 0.1880]];
%set(gcf,'Position', [440, 378, 800, 400]) 
%set(gcf, 'color', 'none');

% plot pore and all zones
for i = 1:numZones
    
    nexttile
    plot(timeSignal, EventSignals(:,i,eventNum)./scale,'Color',colors(i,:),'LineWidth', 2);
    set(gca, 'ylim', [yStart, yEnd], 'FontSize',24);
    set(gca, 'color', 'none');
    
    if i == numZones
        set(gca, 'xlim', [startTime, endTime],'FontSize',24);
    else
        set(gca, 'xlim', [startTime, endTime], 'Xticklabel',[]);
    end
    ax = gca;
    ax.YAxis.Exponent = 0;

end

xlabel(tplot,'Time [s]','FontSize',28);
ylabel(tplot,'Resistance [K\Omega]', 'FontSize',28);
titleTop = sprintf("Event #%d - Resistance\n", eventNum);
titleBottom = sprintf("Cell Type: %s | Device: %s %s | %s", cellType, deviceID, deviceWidth, inputPressure);
titleText = strcat(titleTop, titleBottom);
title(tplot, titleText, 'FontSize',20);

clear tt i titleText1 titleText2 titleText colors tplot ax scale eventNum ...
      timeSignal startTime endTime yStart yEnd

%% *Optional* Plot Frequency Spectrum of Individual Events

eventNum = 5; % choose which event to plot
signal = reshape(EventSignals(:,1:numZones,eventNum),...
                    size(EventSignals(:,1:numZones,eventNum),1),...
                    size(EventSignals(:,1:numZones,eventNum),2));
plotDFT(signal,sampleRate,0,300);

clear eventNum signal

%% Identify Zone Signals from Parsed Events
% Initialize and populate data matrix containing the start and end times of
% each parsed out component (pore, zone 1, zone 2, etc) for each parsed event.

% matrix format:
%{
 each column (9 total): ZonePore    Zone1(in R1)    Zone2(in R2) and so on...   
                        ZonePore_deltaR  ZonePore_baseR  Cell Size
 row 1: start time
 row 2: end time
 each page: different selected event
%}
EventSignalsParsed = zeros(2, numZones+3, size(EventSignals,3));

for eventNum = 1:size(EventSignals,3) % iterate through each selected event
    for z = 1:numZones % for each event, iterate through each zone
        
        % set zone threhold as the median of the resistance signal
        %zt = median(EventSignals(:,z,eventNum));
       
        % set zone threhold as the lower envelope of the resistance signal
        [~,lower] = envelope(EventSignals(:,z,eventNum), sampleRate, 'peak');
        zt = median(lower);
        
    % *** *EDIT* Tunable Parameters *** 
    
        % Threshold value for flagging a sample as an "intersect" via the
        % subtraction method.
        
        % quadratic equation: max value of the event signal in zone 1 as the input.
        % intersectThreshold = (max(EventSignals(:,1,eventNum)))^2*(2e-7);
        
        % alternate equation: n% of the max signal in the respective zone
        %intersectThreshold = zt + ((max(EventSignals(:,z,eventNum))-zt)*0.1);
        intersectThreshold = zt + ((max(rmoutliers(EventSignals(:,z,eventNum)))-zt)*0.05);

        
    % *** Parse Zone Signal Via Baseline Intersects Interval ***
    
        % initialize baseline intersections array 
        intersects = false(length(EventSignals(:,z,eventNum))-1,1);
        
        % subtract (element-wise) the threshold [zt] from the event signal
        sub = abs(EventSignals(:,z,eventNum) - zt);
        
        %{ 
        % find intersections (direct method)
        for j = 1:length(EventSignals(:,z,eventNum))-1
            if (EventSignals(j,z,eventNum) <= zt && ...
                EventSignals(j+1,z,eventNum) >= zt) || ...
               (EventSignals(j,z,eventNum) >= zt && ...
                EventSignals(j+1,z,eventNum) <= zt)
                    intersects(j) = 1;
            end
        end
        %}
        
        % find intersections (subtraction method)
        % * advantage of subtraction method is that what is flagged as an 
        %   "intersect" is tunable via the [intersectThreshold]
        for i = 1:length(sub)
            if sub(i) < intersectThreshold
                intersects(i) = 1;
            end
        end
        intersects = find(intersects==1);
        
        % increase intersect threshold by [factor] if no intersects are found
        factor = 1;
        while isempty(intersects) && factor < 10
            intersects = false(length(EventSignals(:,z,eventNum))-1,1);
            newIntersectThreshold = intersectThreshold + abs(intersectThreshold * factor);
            for j = 1:length(sub)
                if sub(j) < newIntersectThreshold
                    intersects(j) = 1;
                end
            end
            intersects = find(intersects==1);
            factor = factor + 1;
        end
        
        % calculate intersection intervals
        intervals = zeros(length(intersects),1);
        for i = 2:length(intersects)
            intervals(i) = intersects(i) - intersects(i-1);
        end
        
        % the base intersect with the greatest interval is the zone event
        [~, index] = sort(intervals, 'descend');
        if length(index) > 1
            zoneStart = intersects(index(1)-1);
        else
            zoneStart = intersects(index(1));
        end
        zoneEnd = intersects(index(1));
        
    % *** Adjust Start & End Time ***
    
        % #1 - Pore at Beginning of Signal Rule: Pore signal should be in
        % the first half of the entire signal.
        n = 2;
        if z == 1 
            while zoneStart > (length(EventSignals(:,z,eventNum))-1)/2 && ...
                    n < floor(length(index)/2)
                zoneStart = intersects(index(n)-1);
                zoneEnd = intersects(index(n));
                n = n + 1;
            end
        end
    
        % #2 - Start/End Time Linearity Rule: Pore should come before Zone 1,
        %      Zone 1 should come before Zone 2, etc.
        n = 2;
        if z > 2 && z < numZones
            while zoneStart < EventSignalsParsed(1,z-1,eventNum) && ...
                    n < floor(length(index)/2)
                zoneStart = intersects(index(n)-1);
                zoneEnd = intersects(index(n));
                n = n + 1;
            end
        end
    
        % #3 - Real Signal Rule: median signal between zoneStart and zoneEnd should
        %      be greater than the baseline median
        while median(EventSignals(zoneStart:zoneEnd,z,eventNum)) <= zt*1.5 && ...
                n < floor(length(index)/2)
            zoneStart = intersects(index(n)-1);
            zoneEnd = intersects(index(n));
            n = n + 1;
        end
    
        % determine the minimum local min value in the zone
        zoneLocalMin = islocalmin(EventSignals(zoneStart:zoneEnd,z,eventNum),...
            'MinProminence',5000,'MaxNumExtrema',3);
        zoneSignal = (EventSignals(zoneStart:zoneEnd,z,eventNum));
        zoneMinValue = min(zoneSignal(zoneLocalMin==1));
        
        % use the minimum local min value as a new start/end threshold
        if EventSignals(zoneStart,z,eventNum) < zoneMinValue 
            while EventSignals(zoneStart,z,eventNum) < zoneMinValue && z > 1
                zoneStart = zoneStart + 1;
            end
        end
        
        if EventSignals(zoneEnd,z,eventNum) < zoneMinValue 
            while EventSignals(zoneEnd,z,eventNum) < zoneMinValue && z > 1
                zoneEnd = zoneEnd - 1;
            end
        end
        
    % *** Populate [EventSignalsParsed] Matrix ***
     
        EventSignalsParsed(1,z,eventNum) = zoneStart;
        EventSignalsParsed(2,z,eventNum) = zoneEnd;
    end
    
end

clear eventNum z i j intersectThreshold intersects zt sub intervals index ...
      zoneStart zoneEnd zoneLocalMin zoneSignal zoneMinValue n factor ...
      newIntersectThreshold upper

%% Identify and Estimate Pore DeltaR
% Populate [EventSignalParsed] matrix with estimated pore delta R and the
% calculated cell size using the voltage signal from zone 1 (V1). 

% initiate logical for flagging which events to remove when pore is
% misidentified
remove = false(size(EventSignals,3),1);

for eventNum = 1:size(EventSignals,3) % iterate through each selected event
    
% *** Resistance Identification and Cell Size Calculation ***

    % event signal (using Rpore)
    timeSignal = EventSignals(:,numZones*2+1,eventNum);
    poreStartIndex = find(t==timeSignal(1));
    poreEndIndex = find(t==timeSignal(EventSignalsParsed(1,1,eventNum)));
    
    % baseline resistance value
    baseR = median(Rf(poreStartIndex:poreEndIndex,1));
    
    % delta resistance value
    poreR = median(EventSignals(EventSignalsParsed(1,1,eventNum):EventSignalsParsed(2,1,eventNum),1,eventNum));
    poreBase = median(EventSignals(1:EventSignalsParsed(1,1,eventNum),1,eventNum));
    poreDeltaR = poreR-poreBase;
    
    % flag event for removal due to pore misidentifcation when:
    % 1.) pore delta R is negative (poreR is less than baseR)
    if poreDeltaR < 0 
        remove(eventNum) = 1;
    end
    
% *** Populate [EventSignalParsed] Matrix ***
    
    % populate pore deltaR
    EventSignalsParsed(1,numZones+1,eventNum) = poreDeltaR;
    
    % populate baseR
    EventSignalsParsed(1,numZones+2,eventNum) = baseR;
    
    % calculate and store cell diameter [um]
    EventSignalsParsed(1,numZones+3,eventNum) = CellDiameter(poreDeltaR, baseR, Deff, L);
    
end

% filter out events flagged for removal
EventSignals(:,:,remove) = [];
EventSignalsParsed(:,:,remove) = [];

fprintf("Total %d pore identifications unsuccessful.\n", nnz(remove));

clear eventNum remove timeSignal poreStartIndex poreEndIndex ...
      poreR baseR poreDeltaR poreBase

%% *Optional* Plot Individual Events With Parsed Components

figure(1)

tplot = tiledlayout(numZones,1);
tplot.TileSpacing = 'compact';
tplot.Padding = 'compact';
set(gcf,'Position',[50 400 600 600])

% plot parameters
eventNum = 3; % choose which event to plot
timeSignal = EventSignals(:,numZones*2+1,eventNum);
startTime = timeSignal(1); endTime = timeSignal(end);
scale = 10^3;
yStart = min(min(EventSignals(:,1:numZones,eventNum)))/scale;
yEnd = max(max(EventSignals(:,1:numZones,eventNum)))/scale;
colors = [[60/255,60/255,60/255];[0, 0.4470, 0.7410]; ...
          [0.8500, 0.3250, 0.0980] ;[0.9290, 0.6940, 0.1250]; ...
          [0.4940, 0.1840, 0.5560] ;[0.4660, 0.6740, 0.1880]];
zt = zeros(numZones, length(EventSignals(:,:,eventNum)));
for k = 1:numZones
    zt(k,:) = linspace(median(EventSignals(1:1000,k,eventNum)), ...
                       median(EventSignals(end-1000:end,k,eventNum)), ...
                       length(EventSignals(:,k,eventNum)));
end
%set(gcf,'Position', [440, 378, 800, 400]) 
%set(gcf, 'color', 'none');

% plot zones
for i = 1:numZones
    
    nexttile
    
    plot(timeSignal(EventSignalsParsed(1,i,eventNum):EventSignalsParsed(2,i,eventNum)), ...
        EventSignals(EventSignalsParsed(1,i,eventNum):EventSignalsParsed(2,i,eventNum),i,eventNum)./scale,...
        'LineWidth', 4,'Color',colors(i,:)); hold on
    plot(timeSignal, EventSignals(:,i,eventNum)./scale,'Color',[160/255,160/255,160/255],'LineWidth', .5); hold on
    %plot(timeSignal, EventSignals(:,i,eventNum).*0+median(EventSignals(:,i,eventNum)),...
        %'LineWidth', 0.5, 'color', [0.9247, 0.3297, 0.2244]); 
    plot(timeSignal, zt(i,:).'./scale, 'LineWidth', 1, 'Color','red', 'LineStyle', '-');  hold on
    
    % plot pore Delta R
    if i == 0
        poreDeltaR = EventSignalsParsed(1,numZones+1,eventNum);
        plot(timeSignal(EventSignalsParsed(1,i,eventNum):EventSignalsParsed(2,i,eventNum)), ...
        EventSignals(EventSignalsParsed(1,i,eventNum):EventSignalsParsed(2,i,eventNum),i,eventNum).*0+poreDeltaR./scale,...
        'LineWidth', 2,'Color','red')
    end
    
     set(gca, 'ylim', [yStart, yEnd], 'FontSize',16);
    
    if i == numZones
        set(gca, 'xlim', [startTime, endTime],'FontSize',16);
    else
        set(gca, 'xlim', [startTime, endTime], 'Xticklabel',[]);
    end
    ax = gca;
    ax.YAxis.Exponent = 0;

end

xlabel(tplot,'Time [s]','FontSize',20);
ylabel(tplot,'\DeltaR [K\Omega]', 'FontSize',20);
titleTop = sprintf("Event #%d - Resistance\n", eventNum);
titleBottom = sprintf("Cell Type: %s | Device: %s %s | %s", cellType, deviceID, deviceWidth, inputPressure);
titleText = strcat(titleTop, titleBottom);
title(tplot,titleText, 'FontSize',20);

%clear tplot eventNum timeSignal startTime endTime scale yStart yEnd colors ...
     % i titleTop titleBottom titleText ax zt



%% Inititate GUI for Manual QC of Zone Signal & Identification 
% Plot cell transit events for each zone with identified zone signal 

% initiate manual selection GUI
ZoneQCGUI(EventSignals, EventSignalsParsed, numZones, fileDir, fileTestName);

%% Retreive Selected Events & *REMOVE* Non-Selected Events

% update EventSignalsParsed with manual updates
% loads a [nx1] matrix named [EventSignalsParsed] into the workspace
load(strcat(fileDir, fileTestName, '-EventSignalsParsedManual.mat'))

% retrieve selected events 
% loads a [nx1] matrix named [selectedEvents] into the workspace
load(strcat(fileDir, fileTestName, '-EventSelection2.mat'))
fprintf("Total %d total events selected.\n", nnz(selectedEvents2));

% filter out events flagged for removal
EventSignals(:,:,~(selectedEvents2)) = [];
EventSignalsParsed(:,:,~(selectedEvents2)) = [];

%% Calculate TT & Fluid Velocities *SKIP* Quality Control
% Remove events with calculated fluid velocities less than zone velocity

% initialize transit time array
transitTimes = zeros(size(EventSignals,3),numZones);

% calculate transit times
for z = 1:numZones
    for eventNum = 1:length(transitTimes)
        
        % zone start and end sample
        zoneStart = EventSignalsParsed(1,z,eventNum);
        zoneEnd = EventSignalsParsed(2,z,eventNum);
        
        % convert sample to time signal
        time = EventSignals(zoneStart:zoneEnd,numZones*2+1,eventNum);
        
        % transit time (recorded in milliseconds)
        transitTimes(eventNum,z) = (time(end) - time(1)) * 1000 ;
        
    end
end

% calculate cell in contraction channel velocity [m/s]
Vzone = (Lzone ./ transitTimes(:,2:end));

% calculate cell in pore velocity (reference velocity) [m/s] 
Vpore = (Lpore ./ transitTimes(:,1));
    
% calculate fluid velocity in the zones using the reference velocity [m/s] 
Vfluid = Vpore .* (w_np/mean([w_max,w_min]));

% calculate flow veloicty
Vflow = Vfluid - Vzone;

clear z eventNum threshold remove

%% Identify and Estimate Pore DeltaR
% Populate [EventSignalParsed] matrix with estimated pore delta R and the
% calculated cell size using the voltage signal from zone 1 (V1). 

% initiate logical for flagging which events to remove when pore is
% misidentified
remove = false(size(EventSignals,3),1);

for eventNum = 1:size(EventSignals,3) % iterate through each selected event
    
% *** Resistance Identification and Cell Size Calculation ***

    % event signal (using Rpore)
    timeSignal = EventSignals(:,numZones*2+1,eventNum);
    poreStartIndex = find(t==timeSignal(1));
    poreEndIndex = find(t==timeSignal(EventSignalsParsed(1,1,eventNum)));
    
    % baseline resistance value
    baseR = median(Rf(poreStartIndex:poreEndIndex,1));
    
    % delta resistance value
    poreR = median(EventSignals(EventSignalsParsed(1,1,eventNum):EventSignalsParsed(2,1,eventNum),1,eventNum));
    poreBase = median(EventSignals(1:EventSignalsParsed(1,1,eventNum),1,eventNum));
    poreDeltaR = poreR-poreBase;
    
    % flag event for removal due to pore misidentifcation when:
    % 1.) pore delta R is negative (poreR is less than baseR)
    if poreDeltaR < 0 
        remove(eventNum) = 1;
    end
    
% *** Populate [EventSignalParsed] Matrix ***
    
    % populate pore deltaR
    EventSignalsParsed(1,numZones+1,eventNum) = poreDeltaR;
    
    % populate baseR
    EventSignalsParsed(1,numZones+2,eventNum) = baseR;
    
    % calculate and store cell diameter [um]
    EventSignalsParsed(1,numZones+3,eventNum) = CellDiameter(poreDeltaR, baseR, Deff, L);
    
end

% filter out events flagged for removal
EventSignals(:,:,remove) = [];
EventSignalsParsed(:,:,remove) = [];

fprintf("Total %d pore identifications unsuccessful.\n", nnz(remove));

clear eventNum remove timeSignal poreStartIndex poreEndIndex ...
      poreR baseR poreDeltaR poreBase
  
%% Rheological Analysis

% Initialize the data matrix contianing the rheological analysis for each
% parsed zone in each selected cell transit event. 
% matrix format:
%{
 each column (numZones-1 total): Zone1 Zone2 Zone3 etc...
 row 1: calculated transit frequency
 row 2: G'
 row 3: G"
 row 4: prestress
 row 5: friction coeffiecient
 row 6: lsq fit squared 2-norm of the residual
 row 7: transit time (in milliseconds)
 row 8: calculated coefficient of drag
 each page: different parsed cell transit event
%}
Rheology = zeros(8, numZones-1, size(EventSignalsParsed,3));

% debugging arrays
compiledDd = zeros(5000, numZones-1, size(EventSignalsParsed,3));
compiledStress = zeros(5000, numZones-1, size(EventSignalsParsed,3));

% begin rheological analysis for each parsed cell transit event signal
for eventNum = 1:size(EventSignalsParsed,3)
    
    % retrieve calculated cell diameter [um]
    d = EventSignalsParsed(1,numZones+3,eventNum);
    
    % perform rheological analysis for each zone
    for  z = 2:numZones
        
        % zone start and end samples
        % value is relative to the start sample of this specific parsed signal
        zoneStart = EventSignalsParsed(1,z,eventNum);
        zoneEnd = EventSignalsParsed(2,z,eventNum);
        
        % time signal
        time = EventSignals(zoneStart:zoneEnd,numZones*2+1,eventNum);
        
        % samples
        %s = linspace(0,length(time)*sampleRate,length(time));
        s = linspace(0,length(time)/sampleRate,length(time));
        
        % transit frequency [rad/s] (rad/s = 2pi*hz)
        w = (periods(z-1)*2*pi) / (time(end)-time(1));
        
        % contraction channel width [um]
        width = transpose((w_max+w_min)/2 + (w_max-w_min)/2 .* cos(w*s));
        
        % strain cell experiences while in contraction channel
        % E = (d-width) ./ d;
        
        % estimated contraction channel effective diameter [unitless]
        Deff_cont = Deff * sqrt(width./w_np);
    
        % retrieve zone resistance signal [ohms]
        signal = EventSignals(zoneStart:zoneEnd,z,eventNum);
        
        % retrieve zone baseline resistance signal [ohms]
        R = EventSignals(zoneStart:zoneEnd,numZones+z,eventNum);
   
        % calculate diameter of cell in contraction channel [um]
        d_cont= real(CellDiameter(signal, R, Deff_cont, L));
        
    % *** calculate diameter of cell in contact with contraction channel (Dd) [um] ***
        
        Dd = (sqrt( (2./(3*width)) .* (d_cont.^3) ));
        compiledDd(1:zoneEnd-zoneStart+1,z-1, eventNum) = Dd;
        
    % *** calculate stress on cell (StressMu) [Pa] ***
        
        
        % From (JHK, ViscoNPS, iScience, 2019)
        StressMu = ( 2 .* dP(d).*(width).*1e-6 ) ./ ...
                   ( pi.*(Dd.*1e-6) );
        compiledStress(1:zoneEnd-zoneStart+1,z-1, eventNum) = StressMu;
        
    % *** Least Squares Fitting to Determine Rheological Parameters ***
        
        % display progress in command window
        if eventNum == 1 && z == 2
            lineLength = fprintf('Performing Rheological Least Squares Fitting On Event %d of %d...\n',...
                                  eventNum, size(EventSignalsParsed,3));
        elseif z == 2
            fprintf(repmat('\b',1,lineLength))
            lineLength = fprintf('Performing Rheological Least Squares Fitting On Event %d of %d...\n',...
                                  eventNum, size(EventSignalsParsed,3));
        end

        % measured x & y data (curve the model is being fitting to)
        xdata = transpose(s);
        ydata = StressMu;
        
        model = @(x,xdata) ( (x(1) * (1-((w_max+w_min)/(2*d))) + ...
                              x(2) * ((w_max-w_min)/(2*d)) * cos(w*xdata) + ...
                              x(3) * ((w_max-w_min)/(2*d)) * sin(w*xdata)) .* ...
                              x(4) );
        
        % unknown variables initial conditions and search bounds
        x0 = [100, 100, 0, 0.3]; % initial values of [prestress*, G', G", friction]
        lb = [0, 0, 0, 0.3]; % lower search bound
        ub = [inf, inf, inf, 0.5]; % upper search bound
        
        % perform least squares fitting of StressMu to our model:
        [x,resnorm] = lsqcurvefit(model,x0,xdata,ydata,lb,ub,...
                                  optimset('Display','off'));
    
    % *** Populate [Rheology] matrix with calculated unknowns ***
    
        Rheology(1,z-1,eventNum) = periods(z-1) / (time(end)-time(1)); % frequency
        Rheology(2,z-1,eventNum) = x(2); % G'
        Rheology(3,z-1,eventNum) = x(3); % G"
        Rheology(4,z-1,eventNum) = x(1); % young's modulus
        Rheology(5,z-1,eventNum) = x(4); % coefficient of friction
        Rheology(6,z-1,eventNum) = resnorm; % fit quality
                                   % resnorm is the squared 2-norm of the 
                                   % residual at x: sum((fun(x,xdata)-ydata).^2)
        Rheology(7,z-1,eventNum) = (time(end) - time(1)) * 1000; % transit time (recorded as milliseconds)
        Rheology(8,z-1,eventNum) = 0; %coeff_of_drag; % calculated Cd
    end
    
end

fprintf(repmat('\b',1,3));
fprintf('\nRheological Least Squares Fitting Complete.\n');

%clear i z d s w E x R Dd zoneStart zoneEnd time width Deff_cont signal ...
%      d_cont StressMu xdata ydata model x0 lb ub fullDir resorm header data ...
%      outputCSV numPeriods resnorm lineLength

%% Save Rheological Data ***

% create new save folder name
saveFolderName = strcat(fileDir, fileTestName,'-Old-JHKModel-Variable-dP-Analysis');

% if save folder does not exist, create new save folder
if ~exist(saveFolderName, 'dir')
    mkdir(saveFolderName); % create the folder
end

% save the data
save(strcat(saveFolderName, '/Rheology.mat'), 'Rheology');

% *** save rheological data as .CSV file ***

% each row: single cell measurement
% each column: parameter ( G' (zone 1 - 5), G" (zone 1 - 5), Cell Size, 
%                          Transit Time (Zone 1 - 5), Frequencies, etc.)

rz = numZones-1; % number of rheological zones

% initialize header for matrix
header = strings(1, 6*rz+6);

% set headers for output CSV file
for i = 1:rz
    header(i) = sprintf("G1-%d",i);
    header(rz+i) = sprintf("G2-%d",i);
    header(2*rz+i) = sprintf("TransitTime%d",i);
    header(3*rz+1+i) = sprintf("Freq%d",i);
    header(4*rz+1+i) = sprintf("Vflow%d",i);
    header(5*rz+1+i) = sprintf("Cd%d",i);
end

header(2*rz+rz+1) = "TransitTimePore";
header(6*rz+1+1) = "CellType";
header(6*rz+1+2) = "CellSize";
header(6*rz+1+3) = "DeviceID";
header(6*rz+1+4) = "DeviceWidth";
header(6*rz+1+5) = "Test";

% initialize data matrix
data = cell(size(EventSignalsParsed,3), 4*rz+5);

% populate rheological data
for eventNum = 1:size(EventSignalsParsed,3) % for each cell transit event
    
    for z = 1:rz % for each zone in each cell transit event
        
        data(eventNum,z) = num2cell(Rheology(2,z,eventNum)); % storage
        data(eventNum,rz+z) = num2cell(Rheology(3,z,eventNum)); % loss
        data(eventNum,2*rz+z) = num2cell(Rheology(7,z,eventNum)); % transit time
        data(eventNum,3*rz+z+1) = num2cell(Rheology(1,z,eventNum)); % frequency
        data(eventNum,4*rz+z+1) = num2cell(Vflow(eventNum,z)); % flow velocity
        data(eventNum,5*rz+z+1) = num2cell(Rheology(8,z,eventNum)); % Cd
        
    end
    
    data(eventNum,2*rz+rz+1) = num2cell(transitTimes(eventNum,1)); % pore transit time
    data(eventNum,6*rz+1+1) = cellstr(cellType); % cell type
    data(eventNum,6*rz+1+2) = num2cell(EventSignalsParsed(1,numZones+3,eventNum)); % cell size
    data(eventNum,6*rz+1+3) = cellstr(deviceID); % device ID
    data(eventNum,6*rz+1+4) = cellstr(deviceWidth(2:end)); % device width
    data(eventNum,6*rz+1+5) = cellstr(fileTestName(5:end)); % test #
    
end

% compile data with header as cell array
outputCSV = [header; data];

% save cell array as .csv file
writematrix(outputCSV,strcat(saveFolderName, '/RheologicalData.csv'))

fprintf("Analysis Data Saved!\n");

clear rz header i eventNum outputCSV data

%% Clear Workspace

clc; clear;

%%

