function [matchedFrameData, processedPanelOutput] = MatchSerialRecvAndFrameData(serialRecvData, frameData)
    processedPanelOutput = ProcessFlashPanelData(serialRecvData);
    processedPanelOutput.discretizedFramesR.UnityRecvFrame = processedPanelOutput.flipDetailsR.UnityRecvFrame(processedPanelOutput.discretizedFramesR.FlipIndex);
    processedPanelOutput.discretizedFramesL.UnityRecvFrame = processedPanelOutput.flipDetailsL.UnityRecvFrame(processedPanelOutput.discretizedFramesL.FlipIndex);

    
    processedPanelOutput.discretizedFramesR.UnityMatchedFrame = nan(height(processedPanelOutput.discretizedFramesR),1);
    frameData.FlashPanelRStatusMatchingExpected = nan(height(frameData), 1);    
    frameData.DetectedFrameOnset = nan(height(frameData), 1);
    frameData.DiscretizedFrameIndex = nan(height(frameData), 1);
    frameData.DiscretizedFlipIndex = nan(height(frameData), 1);

    
    %Remove the rows that are showing duplicates
    [~, idx, ~] = unique(frameData.Frame, 'stable');
    duplicates = frameData(ismember(1:height(frameData), idx), :);
    frameData = duplicates;

    expectedFrameSequenceIdxsR = FindContinuousFrameSequences(processedPanelOutput.frameDetailsR, 1);
    expectedFrameSequenceIdxsL = FindContinuousFrameSequences(processedPanelOutput.frameDetailsL, 1);

    for iSeq = 1:size(expectedFrameSequenceIdxsR,1)
        [processedPanelOutput.discretizedFramesR] = FindMatchedFrames(processedPanelOutput.discretizedFramesR, processedPanelOutput.flipDetailsR, expectedFrameSequenceIdxsR(iSeq,:), frameData);
    end

    
    matchedFrameData = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput);
    fred = 2;
end

function [discretizedFrames] = FindMatchedFrames(discretizedFrames, flipDetails, sequenceIs, frameData)

    if(sequenceIs(1) == sequenceIs(2))
        discFrameVals = discretizedFrames(sequenceIs(1),:);
    else
        discFrameVals = discretizedFrames(sequenceIs(1):(sequenceIs(2)-1),:);
    end
        flipDetailVals = flipDetails(max(discFrameVals{1,2} - 24, 1) : min(discFrameVals{end,2} + 24, end),:);
    
    
    frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame < flipDetailVals.UnityRecvFrame(end),:);
    
    originalReportedFrameStatus = frameVals.FlashPanelRStatus;
    detectedFrameStatus = discFrameVals.Pattern;
    
    % % Cross-correlation
    % [rs, lags] = xcorr(originalReportedFrameStatus, detectedFrameStatus, 24);
    % 
    % % Find the lag that maximizes the cross-correlation
    % lag = max(lags(rs == max(rs))) + 1;
    % 
    % % Check if there is a need for adjustment
    % if lag > 0
    %     % Adjust reportedFrameStatus based on the lag
    %     reportedFrameStatus = [originalReportedFrameStatus(lag:end); originalReportedFrameStatus(1:lag-1)];
    % else
    %     % Adjust for the case when lag is negative
    %     lag = 0;
    %     reportedFrameStatus = originalReportedFrameStatus;
    % end

    lag = 0;
    reportedFrameStatus = originalReportedFrameStatus;

    % Create a comparison matrix
    try
        comparison = [reportedFrameStatus(1:length(detectedFrameStatus)), detectedFrameStatus];
    catch
        fred = 2;
        comparison = [reportedFrameStatus, detectedFrameStatus(1:length(reportedFrameStatus))];

    end
    % Check for equality
    if ~isequal(comparison(:, 1), comparison(:, 2))
        
        shiftFactor = 0;
        % sequenceIs(2) = sequenceIs(2) - 1;
        breakOccurred = false;  % Initialize the variable
    
        % Keep shifting until the sequences match or lag reaches a limit
        while sum(abs(diff(comparison(1:end-1, :), 1, 2))) ~= 0
            %disp("SHIFTING")
            lag = lag + 1;
    
             if (lag > 24)
                 lag = 1;
                 breakOccurred = true;
                 disp(["LAG " lag]);
                 break;
             end
    
            discFrameVals = discretizedFrames(sequenceIs(1):(sequenceIs(2)-1), :);
            flipDetailVals = flipDetails(max(discFrameVals{1, 2} - 24, 1) : min(discFrameVals{end, 2} + 24, end), :);
    
            frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame <= flipDetailVals.UnityRecvFrame(end), :);
    
            % Update reportedFrameStatus after each shift
            try
            reportedFrameStatus = [originalReportedFrameStatus(lag:end); originalReportedFrameStatus(1:lag-1)];
            catch
                disp("out of range")
            end
            % Update comparison matrix
            comparison = [reportedFrameStatus(1:length(detectedFrameStatus)), detectedFrameStatus];
        end

        if ~breakOccurred
            try
                discretizedFrames.UnityMatchedFrame(sequenceIs(1):(sequenceIs(2)-1)) = frameVals.Frame(lag:lag+length(detectedFrameStatus) - 1);
            catch
                fred = 2;
            end
        end
    
    else
       % disp("EQUALITY")
        discretizedFrames.UnityMatchedFrame(sequenceIs(1):(sequenceIs(2)-1)) = frameVals.Frame(1:length(detectedFrameStatus));
    end
end

function [frameData] = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput)
    % Add a column with row indices to discretizedFramesR
    processedPanelOutput.discretizedFramesR.RowIndices = (1:height(processedPanelOutput.discretizedFramesR))';

    rowsWithMatchingFrames = ~isnan(processedPanelOutput.discretizedFramesR.UnityMatchedFrame);
    discretizedDataWithMatchingFrames = processedPanelOutput.discretizedFramesR(rowsWithMatchingFrames, :);

    for iMatching = 1:height(discretizedDataWithMatchingFrames)
        frameVal = discretizedDataWithMatchingFrames.UnityMatchedFrame(iMatching);
        
        % Find the corresponding row in frameData
        matchingRow = frameData.Frame == frameVal;
        
        % Assign the Correct SynchBox Time to the
        % frameData.DetectedFrameOnset (Using Left Panel)
        frameData.DetectedFrameOnset(matchingRow) = processedPanelOutput.flipDetailsL.CorrectedTime(discretizedDataWithMatchingFrames.FlipIndex(iMatching));
        
        % Assign the original index to the frameData.DiscretizedFrameIndex
        frameData.DiscretizedFrameIndex(matchingRow) = discretizedDataWithMatchingFrames.RowIndices(iMatching);

        % Assign FlipIndex to the frameData.DiscretizedFlipIndex 
        frameData.DiscretizedFlipIndex(matchingRow) = discretizedDataWithMatchingFrames.FlipIndex(iMatching);
        
        % Assign the value of whether or not the flashed pattern matched
        % the expected pattern (Using the Right Panel)
        frameData.FlashPanelRStatusMatchingExpected(matchingRow) = processedPanelOutput.frameDetailsR.Accuracy(iMatching);
    end
end



