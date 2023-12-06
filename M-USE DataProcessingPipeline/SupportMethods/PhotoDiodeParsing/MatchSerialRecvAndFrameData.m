function MatchSerialRecvAndFrameData(serialRecvData, frameData)
processedPanelOutput = ProcessFlashPanelData(serialRecvData);
processedPanelOutput.discretizedFramesR.UnityRecvFrame = processedPanelOutput.flipDetailsR.UnityRecvFrame(processedPanelOutput.discretizedFramesR.FlipIndex);
processedPanelOutput.discretizedFramesR.UnityMatchedFrame = nan(height(processedPanelOutput.discretizedFramesR),1);
frameData.DetectedFrameOnset = nan(height(frameData), 1);
frameData.DiscretizedFrameIndex = nan(height(frameData), 1);
frameData.DiscretizedFlipIndex = nan(height(frameData), 1);

expectedFrameSequenceIdxs = FindContinuousFrameSequences(processedPanelOutput.frameDetailsR, 1);



for iSeq = 1:size(expectedFrameSequenceIdxs,1)
    matchedFrames = FindMatchedFrames(processedPanelOutput, expectedFrameSequenceIdxs(iSeq,:), frameData);
end


fred = 2;


function [frameData, processedPanelOutput] = FindMatchedFrames(processedPanelOutput, sequenceIs, frameData)

discFrameVals = processedPanelOutput.discretizedFramesR(sequenceIs(1):sequenceIs(2),:);
flipDetailVals = processedPanelOutput.flipDetailsR(max(discFrameVals{1,2} - 24, 1) : min(discFrameVals{end,2} + 24, end),:);

frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame < flipDetailVals.UnityRecvFrame(end),:);

reportedFrameStatus = frameVals.FlashPanelRStatus;
detectedFrameStatus = discFrameVals.Pattern;

[rs, lags] = xcorr(reportedFrameStatus, detectedFrameStatus, 24);
lag = lags(rs==max(rs)) + 1;

comparison = [reportedFrameStatus(lag:lag+length(detectedFrameStatus) - 1) detectedFrameStatus];

if ~isequal(comparison(:,1), comparison(:,2))
    if sum(abs(diff(comparison(1:end-1,:),1,2))) == 0
        sequenceIs(2) = sequenceIs(2) - 1;
        
        discFrameVals = processedPanelOutput.discretizedFramesR(sequenceIs(1):sequenceIs(2),:);
        flipDetailVals = processedPanelOutput.flipDetailsR(max(discFrameVals{1,2} - 24, 1) : min(discFrameVals{end,2} + 24, end),:);
        
        frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame < flipDetailVals.UnityRecvFrame(end),:);
        disp("shifted")
    else
        disp("RUH ROH")
    end
else
    disp("EQUALITY")
end

processedPanelOutput.discretizedFramesR.UnityMatchedFrame(sequenceIs(1):sequenceIs(2)) = frameVals.Frame(lag:lag+length(detectedFrameStatus) - 2);
fred = 2;