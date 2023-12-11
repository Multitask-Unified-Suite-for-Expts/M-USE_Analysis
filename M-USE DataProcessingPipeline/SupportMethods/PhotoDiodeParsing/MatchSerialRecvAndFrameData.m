function [matchedFrameData, processedPanelOutput] = MatchSerialRecvAndFrameData(serialRecvData, frameData)
    processedPanelOutput = ProcessFlashPanelData(serialRecvData);
    processedPanelOutput.discretizedFramesR.UnityRecvFrame = processedPanelOutput.flipDetailsR.UnityRecvFrame(processedPanelOutput.discretizedFramesR.FlipIndex);
    processedPanelOutput.discretizedFramesL.UnityRecvFrame = processedPanelOutput.flipDetailsL.UnityRecvFrame(processedPanelOutput.discretizedFramesL.FlipIndex);

    
    processedPanelOutput.discretizedFramesR.UnityMatchedFrame = nan(height(processedPanelOutput.discretizedFramesR),1);
    processedPanelOutput.discretizedFramesL.UnityMatchedFrame = nan(height(processedPanelOutput.discretizedFramesL),1);
    frameData.FlashPanelRValidity = nan(height(frameData), 1);    
    frameData.DetectedFrameOnset = nan(height(frameData), 1);
    frameData.DiscretizedFrameIndex = nan(height(frameData), 1);
    frameData.DiscretizedFlipIndex = nan(height(frameData), 1);

    
    %Remove the rows that are showing duplicates
    [~, idx, ~] = unique(frameData.Frame, 'stable');
    duplicates = frameData(ismember(1:height(frameData), idx), :);
    frameData = duplicates;

    expectedFrameSequenceIdxsR = FindContinuousFrameSequences(processedPanelOutput.frameDetailsR, 1);
    for iSeq = 1:size(expectedFrameSequenceIdxsR,1)
        [processedPanelOutput.discretizedFramesR] = FindMatchedFrames(processedPanelOutput.discretizedFramesR, processedPanelOutput.flipDetailsR, expectedFrameSequenceIdxsR(iSeq,:), frameData);
    end

    [processedPanelOutput.discretizedFramesL] = AssignMatchedFramesToDiscretizedFramesLeft(processedPanelOutput);
    matchedFrameData = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput);
    fred = 2;
end

function [discretizedFramesR] = FindMatchedFrames(discretizedFramesR, flipDetails, sequenceIs, frameData)

    if(sequenceIs(1) == sequenceIs(2))
        discFrameVals = discretizedFramesR(sequenceIs(1),:);
    else
        discFrameVals = discretizedFramesR(sequenceIs(1):(sequenceIs(2)-1),:);
    end
        flipDetailVals = flipDetails(max(discFrameVals{1,2} - 24, 1) : min(discFrameVals{end,2} + 24, end),:);
    
    
    frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame < flipDetailVals.UnityRecvFrame(end),:);
    
    originalReportedFrameStatus = frameVals.FlashPanelRStatus;
    detectedFrameStatus = discFrameVals.Pattern;

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
        while ~isequal(comparison(:, 1), comparison(:, 2))
            %disp("SHIFTING")
            lag = lag + 1;
    
             if (lag > 24)
                 lag = 1;
                 breakOccurred = true;
                 disp(["LAG " lag]);
                 break;
             end
    
            discFrameVals = discretizedFramesR(sequenceIs(1):(sequenceIs(2)-1), :);
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

        
        discretizedFramesR.UnityMatchedFrame(sequenceIs(1):(sequenceIs(2)-1)) = frameVals.Frame(lag:lag+length(detectedFrameStatus) - 1);
            
    else
       % disp("EQUALITY")
        discretizedFramesR.UnityMatchedFrame(sequenceIs(1):(sequenceIs(2)-1)) = frameVals.Frame(1:length(detectedFrameStatus));
    end
end

function [discFramesL] = AssignMatchedFramesToDiscretizedFramesLeft(processedPanelOutput)
    
    % Find the rows in the discretizedFramesR that have a matched Frame 
    validRowsInDiscretizedFramesR = find(~isnan(processedPanelOutput.discretizedFramesR.UnityMatchedFrame));
    
    % Get a reference to the discretizedFramesL
    discFramesL = processedPanelOutput.discretizedFramesL;
    
    % Loop through each valid row in discretizedFramesR
    for iDiscFramesR = 1:height(validRowsInDiscretizedFramesR)
        
        % Find the Flip Index for each row with a matched Frame
        discFramesR = validRowsInDiscretizedFramesR(iDiscFramesR);
        flipIndexR = processedPanelOutput.discretizedFramesR.FlipIndex(discFramesR);
        
        % Find the corresponding Flip Index in flipDetailsL
        flipIndexL = find(processedPanelOutput.flipDetailsL.CorrectedTime == processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR));
        
        % Check if a match was not found
        if(isempty(flipIndexL))
            % If no match is found, skip to the next iteration
            % disp(["Match wasn't found for " processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR) ", moving onto the next frame."])
            % val = [val;processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR)];
            continue;
        end
        
        % Find the index of the corresponding row in discretizedFramesL
        iDiscFramesL = find(processedPanelOutput.discretizedFramesL.FlipIndex == flipIndexL);
        
        % Assign the UnityMatchedFrame value from discretizedFramesR to discretizedFramesL
        discFramesL.UnityMatchedFrame(iDiscFramesL) = processedPanelOutput.discretizedFramesR.UnityMatchedFrame(discFramesR);
    end
end

function [frameData] = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput)
    % Add a column with row indices to discretizedFramesR
    processedPanelOutput.discretizedFramesR.RowIndices = (1:height(processedPanelOutput.discretizedFramesR))';

    rowsWithMatchingFramesR = ~isnan(processedPanelOutput.discretizedFramesR.UnityMatchedFrame);
    discretizedDataRWithMatchingFrames = processedPanelOutput.discretizedFramesR(rowsWithMatchingFramesR, :);

    rowsWithMatchingFramesL = ~isnan(processedPanelOutput.discretizedFramesL.UnityMatchedFrame);
    discretizedDataLWithMatchingFrames = processedPanelOutput.discretizedFramesL(rowsWithMatchingFramesL, :);

    for iDiscRowR = 1:height(discretizedDataRWithMatchingFrames)
        frameVal = discretizedDataRWithMatchingFrames.UnityMatchedFrame(iDiscRowR);
        
        % Find the corresponding row in frameData
        matchingRow = frameData.Frame == frameVal;
        
        % Assign the original index to the frameData.DiscretizedFrameIndex
        frameData.DiscretizedFrameIndex(matchingRow) = discretizedDataRWithMatchingFrames.RowIndices(iDiscRowR);

        % Assign FlipIndex to the frameData.DiscretizedFlipIndex 
        frameData.DiscretizedFlipIndex(matchingRow) = discretizedDataRWithMatchingFrames.FlipIndex(iDiscRowR);
        
        % Assign the value of whether or not the flashed pattern matched
        % the expected pattern (Using the Right Panel)
        frameData.FlashPanelRValidity(matchingRow) = processedPanelOutput.frameDetailsR.Accuracy(iDiscRowR);
    end

    for iDiscRowL = 1:height(discretizedDataLWithMatchingFrames)
        frameVal = discretizedDataLWithMatchingFrames.UnityMatchedFrame(iDiscRowL);
        
        % Find the corresponding row in frameData
        matchingRow = frameData.Frame == frameVal;
        
        % Assign the Correct SynchBox Time to the
        % frameData.DetectedFrameOnset (Using Left Panel)
        frameData.DetectedFrameOnset(matchingRow) = processedPanelOutput.flipDetailsL.CorrectedTime(discretizedDataLWithMatchingFrames.FlipIndex(iDiscRowL));
    end
end



