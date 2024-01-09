function data_b = demodData(data, demodSignal, sampleRate, m, b, a)
% demodulates data at frequency implicit in demodSignal and lowpass filters

% calculate baseband signal
data_d = data .* demodSignal; % demodulate signal
data_dp = [data_d(1 : sampleRate/10,:); data_d; data_d(end-sampleRate/10+1 : end,:)]; % pad signal to reduce filtering artifacts
data_dpf = filtfilt(b, a, data_dp); % lowpass filter padded signal
data_df = data_dpf(sampleRate/10+1 : end-sampleRate/10,:); % unpad lowpass filtered signal
data_dfp = [repmat(data_df(1,:), m*10, 1); data_df; repmat(data_df(end,:), m*10, 1)]; % pad signal to reduce resampling artifacts
data_dfpr = resample(data_dfp, 1, m); % resample padded signal
data_b = data_dfpr(10+1 : end-10, :); % unpad baseband signal

end