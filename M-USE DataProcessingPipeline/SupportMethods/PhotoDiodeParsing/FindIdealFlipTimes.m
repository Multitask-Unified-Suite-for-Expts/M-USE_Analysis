function [correctedFlipTimes, correctedFlipIndices] = FindIdealFlipTimes(t_upsample, upsampledFlipIndices, frameRate, windowSizeSeconds, maxLagSeconds, overlapSeconds, smoothCorrWindowSeconds)

frameDur = 1/frameRate;
upSampleRate = mean(diff(t_upsample));

frameSamples = round(frameDur / upSampleRate);
windowSizeSamples = round(windowSizeSeconds / upSampleRate);
maxLagSamples = round(maxLagSeconds / upSampleRate);
overlapSamples = round(overlapSeconds / upSampleRate);
smoothCorrWindowSamples = round(smoothCorrWindowSeconds / upSampleRate);

rawFlips = zeros(size(t_upsample));
rawFlips(upsampledFlipIndices) = 1;



idealFlips = zeros(size(t_upsample));

idealFlipExactTimes = t_upsample(1):frameDur:t_upsample(end);
idealFlips(FindClosestWithWindow(t_upsample, idealFlipExactTimes)) = 1;


xcorrs = [];
lags = [];
rs = [];
ps = [];
correctedFlips = zeros(size(t_upsample));

for windowStart = 1 : windowSizeSamples - overlapSamples : length(t_upsample)

    windowIndices = windowStart : min(windowStart + windowSizeSamples - 1, length(t_upsample));

    idealSignal = idealFlips(windowIndices);
    rawSignal = rawFlips(windowIndices);
    
    [xc,lag] = xcorr(idealSignal, rawSignal, maxLagSamples);

    smoothCorrs = smooth(xc, smoothCorrWindowSamples);
    %smoothing can screw up due to spikes at the edges, so get rid of those
    trim = (maxLagSamples / 10);
    smoothCorrs([1:trim end-trim+1:end]) = 0;

    iMaxCorr = find(smoothCorrs==max(smoothCorrs), 1);
    bestLag = lag(iMaxCorr);

    xcorrs = [xcorrs; xc(iMaxCorr)];
    lags = [lags; bestLag];

    if bestLag < 0
        adjustedIdeal = [zeros(1, -bestLag) idealSignal(1: end + bestLag)];
        firstFlip = find(adjustedIdeal, 1);
        adjustedIdeal(firstFlip:-frameSamples:1) = 1;

    else
        adjustedIdeal = [idealSignal(bestLag + 1 : end) zeros(1, bestLag)];
        lastFlip = find(adjustedIdeal, 1, 'last');
        adjustedIdeal(lastFlip:frameSamples:end) = 1;
    end

    if windowStart > 1 %merge overlaps of windows
        overlapPrev = correctedFlips(windowStart : windowStart + overlapSamples - 1);
        overlapNew = adjustedIdeal(1:overlapSamples);

        prevFlips = find(overlapPrev);
        newFlips = find(overlapNew);

        if length(prevFlips) ~= length(newFlips)
            disp("RUH ROH")
            fred = 2;
        elseif max(abs(prevFlips - newFlips)) > 30 %3 ms
            %likely there is an issue with an edge flip in one stream but
            %not in the other: so cut off the edges
            trim = 50;

            prevFlips2 = prevFlips(prevFlips > trim & prevFlips < length(overlapPrev) - trim);
            newFlips2 = newFlips(newFlips > trim & newFlips < length(overlapNew) - trim);
            if length(prevFlips2) ~= length(newFlips2)
                disp("RUH ROH")
                fred = 2;
            elseif max(abs(prevFlips2 - newFlips2)) > 30 %3 ms
                disp('YO NOOOOO')
                fred = 2;
            else
                overlapFlips = round(mean([prevFlips2; newFlips2]));
                adjustedIdeal(1:overlapSamples) = 0;
                adjustedIdeal(overlapFlips) = 1;
            end

            %gotta add the edges back in
            prevTrim = prevFlips(prevFlips <= trim | prevFlips >= length(overlapPrev) - trim);
            newTrim = newFlips(newFlips <= trim | newFlips >= length(overlapNew) - trim);
            %make sure they are edges
            if length(prevTrim) > 1
                disp('ARGH');
            elseif length(newTrim) > 1
                disp('AAAARGH');
            elseif abs(newTrim - prevTrim) < 1000
                disp('AAARGHG');
            else
                adjustedIdeal(prevTrim) = 1;
                adjustedIdeal(newTrim) = 1;
            end
        else
      
            overlapFlips = round(mean([prevFlips; newFlips]));
            adjustedIdeal(1:overlapSamples) = 0;
            adjustedIdeal(overlapFlips) = 1;
        end
    end
    correctedFlips(windowIndices) = adjustedIdeal;


end

correctedFlipIndices = find(correctedFlips == 1);
correctedFlipTimes = t_upsample(correctedFlipIndices);

