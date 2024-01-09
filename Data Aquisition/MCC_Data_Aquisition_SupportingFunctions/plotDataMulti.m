function plotDataMulti(src, event)
% Alan Dong
% 2018-03-08

% This callback function plots both raw data and baseband data in real time
% subplot 1: PSD of raw data
% subplot 2: Time domain baseband data
% subplot 3: PSD of baseband data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% get globals
[sampleRate, signalBandwidth, Vpp, ampsPerVolt, demodSignal, nfftRaw, nfftBaseband, inputDuration, windowDuration, lineHandle, figureHandle] = getGlobals;

% get data
t = event.TimeStamps; % get timestamps
AI_032 = event.Data; % get raw signal [V]
% AI0 is V_TIA
% AI3 is V1
% AI2 is V2
V_TIA = AI_032(:,1);
V1 = AI_032(:,2);
V2 = AI_032(:,3);
V3 = AI_032(:,4);
V4 = AI_032(:,5);
V5 = AI_032(:,6);
V6 = AI_032(:,7);
VT = V1 + V2 + V3 + V4 + V5 + V6;
if demodSignal(2)~=demodSignal(1) % cheap test for AC
    V1 = detrend(V1);
    V2 = detrend(V2);
    V3 = detrend(V3);
    V4 = detrend(V4);
    V5 = detrend(V5);
    V6 = detrend(V6);
end
I_TIA = ampsPerVolt * V_TIA;
% add more later (demod stuff)

%%

% set parameters
resampleRate = min([5 * signalBandwidth, sampleRate]); % resample rate [samp/s]
m = sampleRate / resampleRate; % resampling ratio (should be integer)
n = 5; % filter order
[b, a] = butter(n, signalBandwidth / (sampleRate/2)); % design lowpass filter
lambda = 0.2^(inputDuration); % forgetting rate for averaging
% p = 0.1; % padding duration [s]

% calculate baseband signal I_TIA
I_TIA_d = I_TIA .* demodSignal; % demodulate signal (this is now half of Vpp) [V]
I_TIA_dp = [I_TIA_d(1:sampleRate/10); I_TIA_d; I_TIA_d(end-sampleRate/10+1:end)]; % pad signal to reduce filtering artifacts
% I_TIA_dpf = lowpass(I_TIA_dp, signalBandwidth, sampleRate); % lowpass filter padded signal
I_TIA_dpf = filtfilt(b, a, I_TIA_dp); % lowpass filter padded signal
I_TIA_df = I_TIA_dpf(sampleRate/10+1 : end-sampleRate/10); % unpad lowpass filtered signal
I_TIA_dfp = [repmat(I_TIA_df(1), m*10, 1); I_TIA_df; repmat(I_TIA_df(end), m*10, 1)]; % pad signal to reduce resampling artifacts
I_TIA_dfpr = resample(I_TIA_dfp, resampleRate, sampleRate); % resample padded signal
I_TIA_b = I_TIA_dfpr(10+1 : end-10); % unpad baseband signal [V]

% calculate baseband signal V1
V1_d = V1 .* demodSignal; % demodulate signal (this is now half of Vpp) [V]
V1_dp = [V1_d(1:sampleRate/10); V1_d; V1_d(end-sampleRate/10+1:end)]; % pad signal to reduce filtering artifacts
% V1_dpf = lowpass(V1_dp, signalBandwidth, sampleRate); % lowpass filter padded signal
V1_dpf = filtfilt(b, a, V1_dp); % lowpass filter padded signal
V1_df = V1_dpf(sampleRate/10+1 : end-sampleRate/10); % unpad lowpass filtered signal
V1_dfp = [repmat(V1_df(1), m*10, 1); V1_df; repmat(V1_df(end), m*10, 1)]; % pad signal to reduce resampling artifacts
V1_dfpr = resample(V1_dfp, resampleRate, sampleRate); % resample padded signal
V1_b = V1_dfpr(10+1 : end-10); % unpad baseband signal [V]

% calculate baseband signal V2
V2_d = V2 .* demodSignal; % demodulate signal (this is now half of Vpp) [V]
V2_dp = [V2_d(1:sampleRate/10); V2_d; V2_d(end-sampleRate/10+1:end)]; % pad signal to reduce filtering artifacts
% V2_dpf = lowpass(V2_dp, signalBandwidth, sampleRate); % lowpass filter padded signal
V2_dpf = filtfilt(b, a, V2_dp); % lowpass filter padded signal
V2_df = V2_dpf(sampleRate/10+1 : end-sampleRate/10); % unpad lowpass filtered signal
V2_dfp = [repmat(V2_df(1), m*10, 1); V2_df; repmat(V2_df(end), m*10, 1)]; % pad signal to reduce resampling artifacts
V2_dfpr = resample(V2_dfp, resampleRate, sampleRate); % resample padded signal
V2_b = V2_dfpr(10+1 : end-10); % unpad baseband signal [V]

% calculate baseband signal V3
V3_d = V3 .* demodSignal; 
V3_dp = [V3_d(1:sampleRate/10); V3_d; V3_d(end-sampleRate/10+1:end)]; 
V3_dpf = filtfilt(b, a, V3_dp);
V3_df = V3_dpf(sampleRate/10+1 : end-sampleRate/10); 
V3_dfp = [repmat(V3_df(1), m*10, 1); V3_df; repmat(V3_df(end), m*10, 1)]; 
V3_dfpr = resample(V3_dfp, resampleRate, sampleRate); 
V3_b = V3_dfpr(10+1 : end-10); 

% calculate baseband signal V4
V4_d = V4 .* demodSignal; 
V4_dp = [V4_d(1:sampleRate/10); V4_d; V4_d(end-sampleRate/10+1:end)]; 
V4_dpf = filtfilt(b, a, V4_dp);
V4_df = V4_dpf(sampleRate/10+1 : end-sampleRate/10); 
V4_dfp = [repmat(V4_df(1), m*10, 1); V4_df; repmat(V4_df(end), m*10, 1)]; 
V4_dfpr = resample(V4_dfp, resampleRate, sampleRate); 
V4_b = V4_dfpr(10+1 : end-10); 

% calculate baseband signal V5
V5_d = V5 .* demodSignal; 
V5_dp = [V5_d(1:sampleRate/10); V5_d; V5_d(end-sampleRate/10+1:end)]; 
V5_dpf = filtfilt(b, a, V5_dp);
V5_df = V5_dpf(sampleRate/10+1 : end-sampleRate/10); 
V5_dfp = [repmat(V5_df(1), m*10, 1); V5_df; repmat(V5_df(end), m*10, 1)]; 
V5_dfpr = resample(V5_dfp, resampleRate, sampleRate); 
V5_b = V5_dfpr(10+1 : end-10); 

% calculate baseband signal V6
V6_d = V6 .* demodSignal; 
V6_dp = [V6_d(1:sampleRate/10); V6_d; V6_d(end-sampleRate/10+1:end)]; 
V6_dpf = filtfilt(b, a, V6_dp);
V6_df = V6_dpf(sampleRate/10+1 : end-sampleRate/10); 
V6_dfp = [repmat(V6_df(1), m*10, 1); V6_df; repmat(V6_df(end), m*10, 1)]; 
V6_dfpr = resample(V6_dfp, resampleRate, sampleRate); 
V6_b = V6_dfpr(10+1 : end-10); 

% calculate baseband signal impedance magnitudes [Ohm]
Z1 = abs(V1_b) ./ abs(I_TIA_b);
Z2 = abs(V2_b) ./ abs(I_TIA_b);
Z3 = abs(V3_b) ./ abs(I_TIA_b);
Z4 = abs(V4_b) ./ abs(I_TIA_b);
Z5 = abs(V5_b) ./ abs(I_TIA_b);
Z6 = abs(V6_b) ./ abs(I_TIA_b);
Z = Z1 + Z2 + Z3 + Z4 + Z5+ Z6;

% time vector for baseband signal [s]
t_b = t(1:m:end);

%%

% calculate raw PSD
psdRaw_prev = getGlobalPsdRaw; % get previous raw PSD
freqRaw = (0:nfftRaw/2) * sampleRate / nfftRaw; % frequency vector for raw PSD [Hz]
I_TIA_FFT = fft(hanning(length(I_TIA)) .* detrend(I_TIA), nfftRaw); % FFT of windowed raw signal
psdRaw = abs(I_TIA_FFT(1:nfftRaw/2+1)).^2 / sampleRate * nfftRaw; % compute raw PSD
psdRaw(2:end-1) = 2 * psdRaw(2:end-1); % negative frequency correction
psdRaw = (1-lambda) * psdRaw + lambda * psdRaw_prev; % average previous and new raw PSD
setGlobalPsdRaw(psdRaw); % set new raw PSD

% calculate baseband PSD
psdBaseband_prev = getGlobalPsdBaseband; % get previous baseband PSD
freqBaseband = (0:nfftBaseband/2) * resampleRate / nfftBaseband; % frequency vector for baseband PSD [Hz]
Z1_FFT = fft(hanning(length(Z1)) .* detrend(Z1), nfftBaseband); % FFT of windowed baseband signal
psdBaseband = abs(Z1_FFT(1:nfftBaseband/2+1)).^2 / resampleRate * nfftBaseband; % compute baseband PSD
psdBaseband(2:end-1) = 2 * psdBaseband(2:end-1); % negative frequency correction
psdBaseband = (1-lambda) * psdBaseband + lambda * psdBaseband_prev; % average previous and new baseband PSD
setGlobalPsdBaseband(psdBaseband); % set new baseband PSD

%% PLOT
figure(figureHandle); % shift focus to figure

nfft = 1024;
ff = (-nfft/2:nfft/2-1)*sampleRate/nfft;

xend = t_b(end) + t_b(end) - t_b(end-1); % maybe neater plotting?
xstart = xend - windowDuration;

plotting_m = 1;

subplot(3,3,1);
plot(t(1:end), V1(1:end));
xlabel('Time [s]');
ylabel('V1 [V]');
axis tight;
grid on;
grid minor;

subplot(3,3,4);
plot(t(1:end), V2(1:end));
xlabel('Time [s]');
ylabel('V2 [V]');
axis tight;
grid on;
grid minor;

subplot(3,3,7);
plot(t(1:end), V3(1:end));
xlabel('Time [s]');
ylabel('V3 [V]');
axis tight;
grid on;
grid minor;

subplot(3,3,2);
plot(t(1:end), V4(1:end));
xlabel('Time [s]');
ylabel('V4 [V]');
axis tight;
grid on;
grid minor;

subplot(3,3,5);
plot(t(1:end), V5(1:end));
xlabel('Time [s]');
ylabel('V5 [V]');
axis tight;
grid on;
grid minor;

subplot(3,3,8);
plot(t(1:end), V6(1:end));
xlabel('Time [s]');
ylabel('V6 [V]');
axis tight;
grid on;
grid minor;

% plot raw PSD
subplot(3,3,3);
plot(freqRaw/1e3, 10*log10(psdRaw));
% set(gca, 'xscale', 'log');
title('Power Spectral Density Of I_{TIA} Signal');
xlabel('Frequency [kHz]');
ylabel('PSD [A^2/Hz] (dB)');
axis tight;
grid on;
grid minor;
drawnow limitrate;

% plot time domain baseband signal
subplot(3,3,6);
% lineHandle.MaximumNumPoints = windowDuration * resampleRate; % set maximum number of points
% addpoints(lineHandle, t_b, R/1e3); % add baseband signal to plot
plot(t_b, Z1/1e3, t_b, Z2/1e3, t_b, Z3/1e3, t_b, Z4/1e3, t_b, Z5/1e3, t_b, Z6/1e3);

title(sprintf('Baseband Signal (Bandwidth = %d Hz)\nZ = %d kOhm (std = %d kOhm)', ...
      round(signalBandwidth), round(mean(Z)/1e3), round(std(Z)/1e3)));
xlabel('Time [s]', 'FontSize', 12);
ylabel('Impedance Magnitude [kOhm]', 'FontSize', 12);
legend('Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6');
%ylim([4000 6500]);
% set(gca, 'xlim', [xstart xend]); % set x axis to window length
grid on;
grid minor;
drawnow limitrate;

% plot baseband PSD
subplot(3,3,9);
plot(freqBaseband/1e3, 10*log10(psdBaseband));
% set(gca, 'xscale', 'log');
title('Power Spectral Density Of Baseband Resistance Signal');
xlabel('Frequency [kHz]');
ylabel('PSD [Ohm^2/Hz] (dB)');
grid on;
grid minor;
drawnow limitrate;
axis tight;

end