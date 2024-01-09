function setGlobals(sampleRate, signalBandwidth, Vpp, ampsPerVolt, demodSignal, nfftRaw, nfftBaseband, inputDuration, windowDuration, lineHandle, figureHandle)
    global sampleRate_;
    sampleRate_ = sampleRate;
    global signalBandwidth_;
    signalBandwidth_ = signalBandwidth;
    global Vpp_;
    Vpp_ = Vpp;
    global ampsPerVolt_;
    ampsPerVolt_ = ampsPerVolt;
    global demodSignal_;
    demodSignal_ = demodSignal;
    global nfftRaw_;
    nfftRaw_ = nfftRaw;
    global nfftBaseband_;
    nfftBaseband_ = nfftBaseband;
    global inputDuration_;
    inputDuration_ = inputDuration;
    global windowDuration_;
    windowDuration_ = windowDuration;
    global lineHandle_;
    lineHandle_ = lineHandle;
    global figureHandle_;
    figureHandle_ = figureHandle;
end