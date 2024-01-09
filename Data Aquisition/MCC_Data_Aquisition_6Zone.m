%% NPS acquisition script for MCC DAQ
% Alan Dong + Edits (Andre Lai)
% 2021-11-02

% This script begins outputing an analog waveform and collecting data at the
% specified sampling frequency when it is run. When the stop button is clicked,
% it will terminate the session. It does not save the data by default.

%% stop and release previous session if needed
if exist('s', 'var')
    stop(s);
    release(s);
end
fclose all;
close all;
clear;

%% Specify *NEW* Save Folder

saveFolderName = 'example data - example test'; % ***specify the new save folder name

if ~exist(saveFolderName, 'dir')
    mkdir(saveFolderName); % create the folder
end

%% ACQUISITION PARAMETERS
addpath('MCC_Data_Aquisition_SupportingFunctions');
sampleRate = 10e3; % sampling rate [samp/s] // MAKE THIS A MULTIPLE OF 10!
signalBandwidth = 500; % signal bandwidth [Hz] (this sets the lowpass cutoff) // MAKE THIS A SIMPLE FRACTION OF THE SAMPLING RATE!
nfftRaw = min(2^nextpow2(sampleRate/2), 2^15); % FFT size for raw PSD estimate
nfftBaseband = min(2^nextpow2(5*signalBandwidth), 2^15); % FFT size for baseband PSD estimate
inputDuration = 1.0; % length of each input data read [s] (also sets plot update rate) // MAKE THIS NO LESS THAN 0.1
windowDuration = 1.0; % length of plot window [s]
ampsPerVolt = 1e-7; % current preamplifier sensitivity [A/V]

% generate output signal
freqVin = 0; % frequency of sine excitation [Hz] // MAKE THIS NO MORE THAN 1/5th OF THE SAMPLING RATE AND A MULTIPLE OF 10!
%freqVin = 0; % FOR DC ACQUISITION
Vpp = 5; % peak to peak voltage of excitation signal [V]
if freqVin == 0, Vpp = Vpp*2; end % amplitude fix for DC
outputDuration = 1; % length of each output data write [s] (1 second is fine)
outputData = Vpp/2 * cos(2*pi*freqVin*(0:outputDuration*sampleRate-1)/sampleRate).'; % output data to write to output queue
demodSignal = 2 * exp(-1i*2*pi*freqVin*(0:inputDuration*sampleRate-1)/sampleRate).'; % signal same length as the input data read for demodulation

% initialize PSD estimates
psdRaw = zeros(nfftRaw/2+1, 1);
psdBaseband = zeros(nfftBaseband/2+1, 1);

% fix resistance calculation...
%if freqVin == 0, Vpp = Vpp*2; end

% initialize figure
figureHandle = figure; % open a new figure to plot in
set(figureHandle, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.9 0.9]); % almost fullscreen
subplot(3,2,1);
subplot(3,2,2);
subplot(3,2,3);
subplot(3,2,4);
subplot(3,2,5);
subplot(3,2,6);
lineHandle = animatedline;

% set global variables so plotData callback function can access them
setGlobalPsdRaw(psdRaw);
setGlobalPsdBaseband(psdBaseband);
setGlobals(sampleRate, ...
           signalBandwidth, ...
           Vpp, ...
           ampsPerVolt, ...
           demodSignal, ...
           nfftRaw, ...
           nfftBaseband, ...
           inputDuration, ...
           windowDuration, ...
           lineHandle, ...
           figureHandle);

%% ACQUISITION INITIALIZATION

% create and configure data acquisition session
s = daq.createSession('mcc');
addAnalogOutputChannel(s, 'Board0', 0, 'Voltage');
% addAnalogInputChannel(s, 'Board0', [0,3,2], 'Voltage');
addAnalogInputChannel(s, 'Board0', [4,12,13,6,7,0,1], 'Voltage');

% set session parameters
s.Rate = sampleRate; % sampling rate [samp/s]
s.IsContinuous = true;
s.NotifyWhenDataAvailableExceeds = inputDuration * sampleRate; % this sets the input data buffer read rate
s.NotifyWhenScansQueuedBelow = 5 * outputDuration * sampleRate; % this sets when to write to the output queue

% open binary data log file
fid1 = fopen('log.bin', 'w');

% add callback function listeners
btn = uicontrol('Style', 'pushbutton', 'String', 'Stop', 'Position', [10 10 40 20], ...
                'Callback', 'stop(s); release(s)'); % this push button is supposed to stop the analog input object
lh1 = addlistener(s, 'DataAvailable', @(src, event) plotDataMulti(src, event)); % plots data for every input data buffer read
lh2 = addlistener(s, 'DataAvailable', @(src, event) logDataMulti(src, event, fid1)); % logs data to file for every input data buffer read
lh3 = addlistener(s, 'DataRequired', @(src, event) src.queueOutputData(outputData)); % writes output data whenever output queue is low

%% START ACQUISITION

queueOutputData(s, repmat(outputData, round(10/outputDuration), 1)); % start with 10 seconds of output data
fprintf('Output samples in queue = %d\n', s.ScansQueued);
prepare(s);
pause(1); % just in case

s.startBackground(); % this starts the continuous acquisition
fprintf('Acquisition started!\n');
fprintf('Output samples in queue = %d\n', s.ScansQueued);

% loop for printouts
tic;
while s.IsRunning
    pause(5);
    fprintf('Elapsed time = %d seconds\n', round(toc));
    fprintf('Samples acquired = %d\n', s.ScansAcquired);
    fprintf('Output samples in queue = %d\n', s.ScansQueued);
end

% after completion, delete listeners
delete(lh1);
delete(lh2);
delete(lh3);
fclose(fid1);
fprintf('Acquisition stopped!\n');

%% Load data into MATLAB from binary log

fid2 = fopen('log.bin', 'r');
[data, count] = fread(fid2, inf, 'double'); % load data into MATLAB
fclose(fid2);

delete log.bin

%% Save the data to the specified save folder in a .mat file

data = data.';
filename = input('Input file name:', 's');
saveFile = strcat(saveFolderName, '/', filename, '.mat');
i = 1;

while isfile(saveFile)
    saveFile = strcat(saveFolderName, '/', filename, '-',num2str(i), '.mat');
    i = i + 1;
end

save(saveFile, 'data', 'demodSignal', 'ampsPerVolt', 'sampleRate', 'Vpp', 'freqVin');

return
