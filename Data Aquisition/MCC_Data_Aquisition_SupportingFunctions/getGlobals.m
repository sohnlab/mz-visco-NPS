function [sampleRate, signalBandwidth, Vpp, ampsPerVolt, demodSignal, nfftRaw, nfftBaseband, inputDuration, windowDuration, lineHandle, figureHandle] = getGlobals()
    global sampleRate_;
    sampleRate = sampleRate_;
    global signalBandwidth_;
    signalBandwidth = signalBandwidth_;
    global Vpp_;
    Vpp = Vpp_;
    global ampsPerVolt_;
    ampsPerVolt = ampsPerVolt_;
    global demodSignal_;
    demodSignal = demodSignal_;
    global nfftRaw_;
    nfftRaw = nfftRaw_;
    global nfftBaseband_;
    nfftBaseband = nfftBaseband_;
    global inputDuration_;
    inputDuration = inputDuration_;
    global windowDuration_;
    windowDuration = windowDuration_;
    global lineHandle_;
    lineHandle = lineHandle_;
    global figureHandle_;
    figureHandle = figureHandle_;
end