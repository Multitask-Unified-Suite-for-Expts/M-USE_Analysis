function [troughToPeakTransitions, peakToTroughTransitions, ...
    allTransitions] = FindFrameTransitions(lightData, ...
    windowSize, stepSize)

%%
%{
Detects the start and end of frames in data collected from a light sensor
over a patch assumed to be flickering from black to white. It is agnostic
about the flicker rate and the data collection rate.

lightData: Matrix of data from the sensor. Assumes it is arranged
vertically, which each 

windowSize: size (in data points) of the sliding window to use for peak
calculation

troughToPeakTransitions - indices of the start of a rising edge from a 
black frame to a white frame (assuming higher numbers indicate brightness)

peakToTroughTransitions - indices of the start of a falling edge from a 
white to a black frame

allTransitions: the sorted union of troughToPeakTransitions and
peakToTroughTransitions
%}

%% Find peaks

%make lightdata nx1 (assumed for some computations)
if height(lightData) < width(lightData)
    lightData = lightData';
end

if width(lightData) ~= 1
    error('FindFrameTransitions expects lightData to be in vector form.')
end

troughToPeakTransitions = [];
peakToTroughTransitions = [];
for windowStart = 1:stepSize-1:length(lightData)

    windowData = lightData(windowStart : min(windowStart + windowSize-1, end), 1);
    [peaks, iPeaks] = findpeaks(windowData);
    [troughs, iTroughs] = findpeaks(-windowData);

    %check for plateaus
    for i = 1:length(peaks)
        peak = iPeaks(i);
        if (windowData(peak + 1) == windowData(peak))
            plateauEnd = find(windowData(peak:end) ~= windowData(peak), 1) + peak - 1;
            iPeaks(i) = plateauEnd - 1;
        end
    end
    %check for plateaus
    for i = 1:length(troughs)
        peak = iTroughs(i);
        if (windowData(peak + 1) == windowData(peak))
            plateauEnd = find(windowData(peak:end) ~= windowData(peak), 1) + peak - 1;
            iTroughs(i) = plateauEnd - 1;
        end
    end

    %find local threshold,
    threshold = mean(windowData);
    % windowMeans = [windowMeans; windowStart + round(height(windowData)/2) threshold]; %#ok<AGROW>

    %combine peaks and troughs, sort them by time point (index)
    combined = sortrows([iPeaks, peaks; iTroughs, windowData(iTroughs)]);

    %apply rule: trough to peak is a trough below threshold followed by a 
    % peak above threshold, peak to trough is the reverse
    ttop = combined(combined(1:end-1,2) < threshold & combined(2:end,2) > threshold,:);
    ptot = combined(combined(1:end-1,2) > threshold & combined(2:end,2) < threshold,:);

    troughToPeakTransitions = [troughToPeakTransitions; ttop(:,1) + windowStart - 1, windowData(ttop(:,1),:), threshold * ones(height(ttop),1)]; %#ok<AGROW>
    peakToTroughTransitions = [peakToTroughTransitions; ptot(:,1) + windowStart - 1, windowData(ptot(:,1),:), threshold * ones(height(ptot),1)]; %#ok<AGROW>
    
end

%eliminate duplicates
[~,uniqueTtoP,~] = unique(troughToPeakTransitions(:,1));
[~,uniquePtoT,~] = unique(peakToTroughTransitions(:,1));
troughToPeakTransitions = troughToPeakTransitions(uniqueTtoP,:);
peakToTroughTransitions = peakToTroughTransitions(uniquePtoT,:);


allTransitions = sortrows([troughToPeakTransitions; peakToTroughTransitions]);


