function [frameDetailsL, exceptionDetailsL, frameDetailsR, exceptionDetailsR] = ProcessFlashPanelData(serialRecvData)
%%
%{

Takes the analog data recorded from light sensors placed over two flashing 
panels, and does various clever timing checks with them

analogData: the Analog field of an M-USE SerialRecvData folder, 

%}

%% Constants

frameRate = 60;
expectedFrameDur = 1/frameRate;

upsampleRate = 10000;

expectedLeft = [0 1]';
expectedRight = [0 0 0  0 0 1  0 1 0  0 1 1  1 0 0  1 0 1  1 1 0  1 1 1]';

%% Find Frame Onsets, Flip to Flip Details


disp("Finding frame onsets in raw data")


analog_orig = serialRecvData.Analog(:,["UnityRecvFrame", "SynchBoxTime", "LightL", "LightR"]);
disp("    Left signal")
[flipDetailsL, flipSummaryL] = GetFrameAnalogDetails(analog_orig, "LightL", expectedFrameDur, 600);
disp(flipSummaryL);
disp("    Right signal")
[flipDetailsR, flipSummaryR] = GetFrameAnalogDetails(analog_orig, "LightR", expectedFrameDur, 600);
disp(flipSummaryR);

%% Find Closest Upsampled Frame Onsets, Corrected Flip To Flip Details


disp("Finding upsampled flip times")
t_orig = serialRecvData.Analog.SynchBoxTime;
t_upsample = t_orig(1) : 1/upsampleRate : t_orig(end);
upsampledFlipIndices = FindClosestWithWindow(t_upsample, flipDetailsL.SynchBoxTime);

[correctedFlipTimes, ~] = FindIdealFlipTimes(t_upsample, upsampledFlipIndices, frameRate, 20, 0.02, 1, 0.002);

closestIndicesL = FindClosestWithWindow(correctedFlipTimes, flipDetailsL.SynchBoxTime);
closestIndicesR = FindClosestWithWindow(correctedFlipTimes, flipDetailsR.SynchBoxTime);

flipDetailsL.CorrectedTime = correctedFlipTimes(closestIndicesL)';
flipDetailsR.CorrectedTime = correctedFlipTimes(closestIndicesR)';

%find out proportions of detected flips that aren't close to corrected
%times

flipDetailsL.NumFrames = [diff(closestIndicesL); NaN];
flipDetailsR.NumFrames = [diff(closestIndicesR); NaN];

discretizedFramesL = GenerateDiscretizedFrameVector(flipDetailsL);
discretizedFramesR = GenerateDiscretizedFrameVector(flipDetailsR);


[frameDetailsL, exceptionDetailsL] = FlashPanelPatternProcessor(discretizedFramesL, expectedLeft);
[frameDetailsR, exceptionDetailsR] = FlashPanelPatternProcessor(discretizedFramesR, expectedRight);
% problemFramesL = DetectPatternViolations(discretizedFramesL, expectedLeft, 600);
% problemFramesR = DetectPatternViolations(discretizedFramesR, expectedRight, 600);


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
analogFrameData.Threshold = allTransitions(:,3);
analogFrameData.FlashPanelStatus = analogFrameData.(column) < allTransitions(:,3);
analogFrameData.OneMsOff = [abs(expectedFrameDur - duration) > 0.001; NaN];
analogFrameData.FiveMsOff = [abs(expectedFrameDur - duration) > 0.005; NaN];


frameSummary.MeanFtoFDuration = mean(duration);
frameSummary.MedianFtoFDuration = median(duration);
frameSummary.FtoFStdDev = std(duration);
frameSummary.FtoFStdErr = std(duration) ./ (length(duration) - 1) ^ 2;
frameSummary.NumFlips = length(duration);
frameSummary.NumGoodFlips = length(duration) - sum(analogFrameData.OneMsOff, 'omitnan');
frameSummary.NumFlipsOneMsOff = sum(analogFrameData.OneMsOff, 'omitnan');
frameSummary.NumFlipsFiveMsOff = sum(analogFrameData.FiveMsOff, 'omitnan');
frameSummary.PropGoodFlips = frameSummary.NumGoodFlips / length(duration);
frameSummary.PropFlipsOneMsOff = frameSummary.NumFlipsOneMsOff / length(duration);
frameSummary.PropFlipFiveMsOff = frameSummary.NumFlipsFiveMsOff / length(duration);




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



function FrameOnsets = findIdealFrameOnsets(analogData, frameRate, sampleRate, windowSizeSeconds, maxLag)
%%
%{
frameRate - of the monitor
sampleRate - of the light sensor data
windowSizeFrames - to compute xcorrs over, in frames (not s)
maxLag - max lag for xcorr in s
%}

frameDur = 1 / frameRate;


windowSizeFrames = windowSizeSeconds / frameDur;

idealTemplate = zeros(ceil(windowSizeSeconds * 1 / sampleRate), 1);

windowSizeSamples = windowSizeSeconds * sampleRate;
frameSizeSamples = round(frameDur * sampleRate);


bToW = [zeros(1, frameSizeSamples) ones(1, frameSizeSamples)];

idealTemplate = repmat(bToW, 1, ceil(windowSizeFrames/2));

t = analogData.SynchBoxTime;
% tnorm = t - t(1);
% end

rs = [];
lags = [];
matches = [];
correctedStatus = [];
maxLagSamples = maxLag / (1 / sampleRate);


for timer = t(1) : windowSizeSeconds : t(end)
    windowStart = find(t > timer, 1) - 1;
    measuredWindow = analogData.FlashPanelStatus(windowStart : min(windowStart + length(idealTemplate) - 1, end));
    windowT = analogData.SynchBoxTime(windowStart : min(windowStart + length(idealTemplate) - 1, end));
    windowTNorm = windowT - windowT(1);

    %generate accurately-timed ideal sample
    timer = windowTNorm(1);

    idealOnsets = windowTNorm(1):frameDur:windowTNorm(end);
    idealSignal = zeros(length(measuredWindow), 1);
    for iFrame = 1:2:windowSizeFrames
        idealSignal(windowTNorm >= idealOnsets(iFrame) & windowTNorm < idealOnsets(iFrame + 1)) = 1;
    end

    [r, lag] = xcorr(idealSignal, measuredWindow, maxLagSamples);
    % display(max(r))
    matches = [matches; min(lag(r==max(r)))];
end

meanMatch = round(mean(matches));


fred = 2;