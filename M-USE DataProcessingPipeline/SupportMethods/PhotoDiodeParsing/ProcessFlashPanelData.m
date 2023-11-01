function ProcessFlashPanelData(serialRecvData)
%%
%{

Takes the analog data recorded from light sensors placed over two flashing 
panels, and does various clever timing checks with them

analogData: the Analog field of an M-USE SerialRecvData folder, 

%}
%%

frameRate = 60;
expectedFrameDur = 1/frameRate;

upsampleRate = 10000;

t_orig = serialRecvData.Analog.SynchBoxTime;
t_upsample = t_orig(1) : 1/upsampleRate : t_orig(end);

analog_orig = serialRecvData.Analog(:,["UnityRecvFrame", "SynchBoxTime", "LightL", "LightR"]);

analog_upsample = array2table(interp1(t_orig, table2array(analog_orig(:,2:4)), t_upsample), 'VariableNames', ["SynchBoxTime", "LightL", "LightR"]);

[analogFrameDetailsL, frameSummaryL] = GetFrameAnalogDetails(analog_orig, "LightL", expectedFrameDur, 600);
[analogFrameDetailsR, frameSummaryR] = GetFrameAnalogDetails(analog_orig, "LightR", expectedFrameDur, 600);


[analogFrameDetailsL_upsample, frameSummaryL_upsample] = GetFrameAnalogDetails(analog_upsample, 'LightL', expectedFrameDur, 20000);

analog_upsample.FlashPanelStatus = GetUpsampledFlashPanelStatus(height(analog_upsample), analogFrameDetailsL_upsample);

FrameOnsets = findIdealFrameOnsets(analog_upsample, 60, 10000, 100, 0.02);


fred = 2;



function [analogFrameData, frameSummary] = GetFrameAnalogDetails(sData, column, expectedFrameDur, windowSize)

%windowSize = 600; %duration of thresholding window in samples
stepSize = 500; 

lightData = sData.(column);
[~, ~, allTransitions] = FindFrameTransitions(lightData, windowSize, stepSize);

analogFrameData = sData(allTransitions(:,1),:);
analogFrameData.AnalogIndex  = allTransitions(:,1);

time = analogFrameData.SynchBoxTime;
duration = diff(time);
analogFrameData.Duration = [duration; NaN];

frameSummary.MeanFrameDuration = mean(duration);
frameSummary.MedianFrameDuration =median(duration);

analogFrameData.Threshold = allTransitions(:,3);
analogFrameData.FlashPanelStatus = analogFrameData.(column) < allTransitions(:,3);
analogFrameData.ProblemFrame = [abs(expectedFrameDur - duration) > 0.5; NaN];



function logicalLight = GetUpsampledFlashPanelStatus(n_upsample, frameDetails)

logicalLight = nan(n_upsample,1);

for iFrame = 1:height(frameDetails)-1
    logicalLight(frameDetails.AnalogIndex(iFrame) : frameDetails.AnalogIndex(iFrame + 1) - 1) = frameDetails.FlashPanelStatus(iFrame);
end

if frameDetails.AnalogIndex(1) ~= 1
    logicalLight(1:frameDetails.AnalogIndex(1)) = abs(frameDetails.FlashPanelStatus(1) - 1);
end

if frameDetails.AnalogIndex(end) ~= height(frameDetails)
    logicalLight(frameDetails.AnalogIndex(end) : end) = abs(frameDetails.FlashPanelStatus(end));
end



function FrameOnsets = findIdealFrameOnsets(analogData, frameRate, sampleRate, windowSizeFrames, maxLag)
%%
%{
frameRate - of the monitor
sampleRate - of the light sensor data
windowSizeFrames - to compute xcorrs over, in frames (not s)
maxLag - max lag for xcorr in s
%}

frameDur = 1 / frameRate;

t = analogData.SynchBoxTime;
tnorm = t - t(1);

windowSizeSeconds = windowSizeFrames * frameDur;

idealTemplate = zeros(ceil(windowSizeSeconds * 1 / sampleRate), 1);
for timer = 0 : frameDur * 2 : windowSizeSeconds - frameDur * 2

    idealTemplate(tnorm >= timer & tnorm < timer + frameDur) = 1;

end

rs = [];
lags = [];
matches = [];
maxLagSamples = maxLag / (1 / sampleRate);
for timer = 0 : windowSizeSeconds : tnorm(end)
    windowStart = find(tnorm > timer, 1) - 1;
    measuredWindow = analogData.FlashPanelStatus(windowStart : min(windowStart + length(idealTemplate) - 1, end));
    [r, lag] = xcorr(idealTemplate, measuredWindow, maxLagSamples);
    % display(max(r))
    matches = [matches; min(lag(r==max(r)))];
end

meanMatch = round(mean(matches));


fred = 2;