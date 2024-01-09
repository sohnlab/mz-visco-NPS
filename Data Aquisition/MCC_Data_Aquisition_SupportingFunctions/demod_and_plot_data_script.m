% load('test5.mat');
%load('test9.mat');

% set parameters
signalBandwidth = 50;
resampleRate = min([5 * signalBandwidth, sampleRate]); % resample rate [samp/s]
m = sampleRate / resampleRate; % resampling ratio (should be integer)
n = 5; % filter order
[b, a] = butter(n, signalBandwidth / (sampleRate/2)); % design lowpass filter
% p = 0.1; % padding duration [s]

% calculate baseband signal I_TIA
data = data.';
I_TIA = data(1:3:end) / 10e6;
V1 = data(2:3:end);
V2 = data(3:3:end);

demodSignal = repmat(demodSignal, length(I_TIA)/length(demodSignal), 1);

I_TIA_d = I_TIA .* demodSignal; % demodulate signal (this is now half of Vpp) [V]
I_TIA_dp = [I_TIA_d(1:sampleRate/10); I_TIA_d; I_TIA_d(end-sampleRate/10+1:end)]; % pad signal to reduce filtering artifacts
% I_TIA_dpf = lowpass(I_TIA_dp, signalBandwidth, sampleRate); % lowpass filter padded signal
I_TIA_dpf = filtfilt(b, a, I_TIA_dp); % lowpass filter padded signal
I_TIA_df = I_TIA_dpf(sampleRate/10+1 : end-sampleRate/10); % unpad lowpass filtered signal
I_TIA_dfp = [repmat(I_TIA_df(1), m*10, 1); I_TIA_df; repmat(I_TIA_df(end), m*10, 1)]; % pad signal to reduce resampling artifacts
I_TIA_dfpr = resample(I_TIA_dfp, resampleRate, sampleRate); % resample padded signal
I_TIA_b = I_TIA_dfpr(10+1 : end-10);

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

% calculate baseband signal impedance magnitudes [Ohm]
Z1 = abs(V1_b) ./ abs(I_TIA_b);
Z2 = abs(V2_b) ./ abs(I_TIA_b);
Z = Z1 + Z2;

figure(1); clf; hold on;
plot(Z1);
plot(Z2);


