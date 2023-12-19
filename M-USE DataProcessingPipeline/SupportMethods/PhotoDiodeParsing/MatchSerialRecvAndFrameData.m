function [matchedFrameData, processedPanelOutput] = MatchSerialRecvAndFrameData(serialRecvData, frameData)
    processedPanelOutput = ProcessFlashPanelData(serialRecvData);
    processedPanelOutput.discretizedFramesR.UnityRecvFrame = processedPanelOutput.flipDetailsR.UnityRecvFrame(processedPanelOutput.discretizedFramesR.FlipIndex);
    processedPanelOutput.discretizedFramesL.UnityRecvFrame = processedPanelOutput.flipDetailsL.UnityRecvFrame(processedPanelOutput.discretizedFramesL.FlipIndex);

    
    processedPanelOutput.discretizedFramesR.UnityMatchedFrame = nan(height(processedPanelOutput.discretizedFramesR),1);
    processedPanelOutput.discretizedFramesL.UnityMatchedFrame = nan(height(processedPanelOutput.discretizedFramesL),1);
    frameData.FrameValidity = nan(height(frameData), 1);    
    frameData.FrameOnsetSyncBoxTime = nan(height(frameData), 1);
    frameData.DiscretizedFrameIndex = nan(height(frameData), 1);
    frameData.FlipIndex = nan(height(frameData), 1);

    
    %Remove the rows that are showing duplicates
    frameData(diff(frameData.Frame) == 0, :) = [];
    % [~, idx, ~] = unique(frameData.Frame, 'stable');
    % nonduplicates = frameData(ismember(1:height(frameData), idx), :);
    % frameData = nonduplicates;

    validFrameSequenceIdxsR = FindValidFrameSequences(processedPanelOutput.frameDetailsR, 1);
    for iSeq = 1:size(validFrameSequenceIdxsR,1)
        [processedPanelOutput.discretizedFramesR] = FindMatchedFrames(processedPanelOutput.discretizedFramesR, processedPanelOutput.flipDetailsR, validFrameSequenceIdxsR(iSeq,:), frameData);
    end

    [processedPanelOutput.discretizedFramesL] = AssignMatchedFramesToDiscretizedFramesLeft(processedPanelOutput);
    matchedFrameData = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput);
    fred = 2;
end

function [detectedFrames] = FindMatchedFrames(detectedFrames, detectedFlipDetails, detectedFrameBoundaries, unityReportedFrameData)


    %find the frames in frameData that corresponed to the sequence in
    %DiscretizedFrames that correspond

    %exact window of detected & discretzed frames
    validDetectedFrameWindow = detectedFrames(detectedFrameBoundaries(1):detectedFrameBoundaries(2),:);

    %flipDetails contains frames in which unity received syncbox data
    flipDetailVals = detectedFlipDetails(max(validDetectedFrameWindow{1,2} - 24, 1) : min(validDetectedFrameWindow{end,2} + 24, end),:);
    frameVals = unityReportedFrameData(unityReportedFrameData.Frame <= flipDetailVals.UnityRecvFrame(1) & unityReportedFrameData.Frame < flipDetailVals.UnityRecvFrame(end),:);

    originalReportedFrameStatus = frameVals.FlashPanelRStatus;
    detectedFrameStatus = validDetectedFrameWindow.Status;

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
                 breakOccurred = true;
                 disp(["LAG " lag]);
                 lag = 1;

                 break;
             end
    
            % discFrameVals = discretizedFramesR(sequenceIs(1):(sequenceIs(2)-1), :);
            % flipDetailVals = flipDetails(max(discFrameVals{1, 2}, 1) : min(discFrameVals{end, 2} + 24, end), :);
            % 
            % frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame <= flipDetailVals.UnityRecvFrame(end), :);
    
            % Update reportedFrameStatus after each shift
            try
                reportedFrameStatus = [originalReportedFrameStatus(lag:end); originalReportedFrameStatus(1:lag-1)];
            catch
                disp("out of range")
            end
            % Update comparison matrix
            comparison = [reportedFrameStatus(1:length(detectedFrameStatus)), detectedFrameStatus];
        end

        try
            detectedFrames.UnityMatchedFrame(detectedFrameBoundaries(1):(detectedFrameBoundaries(2)-1)) = frameVals.Frame(lag:lag+length(detectedFrameStatus) - 1);
        catch
            detectedFrames.UnityMatchedFrame(detectedFrameBoundaries(1):detectedFrameBoundaries(1)+height(frameVals.Frame(lag:end))-1) = frameVals.Frame(lag:end);

            fred = 2;
        end
    else
       % disp("EQUALITY")
       % Determine the valid range based on the lengths of sequence and frameVals.Frame
    range = detectedFrameBoundaries(1):(detectedFrameBoundaries(1) + length(frameVals.Frame(1:length(detectedFrameStatus))) - 1);

    % Update only the existing rows in discretizedFramesR.UnityMatchedFrame
    detectedFrames.UnityMatchedFrame(range) = frameVals.Frame(1:length(detectedFrameStatus));

    end
end

function [discFramesL] = AssignMatchedFramesToDiscretizedFramesLeft(processedPanelOutput)
    
    % Find the rows in the discretizedFramesR that have a matched Frame 
    validRowsInDiscretizedFramesR = processedPanelOutput.discretizedFramesR(~isnan(processedPanelOutput.discretizedFramesR.UnityMatchedFrame), :);
    
    % Get a reference to the discretizedFramesL
    discFramesL = processedPanelOutput.discretizedFramesL;
    
    % Loop through each valid row in discretizedFramesR
    for iDiscFramesR = 1:height(validRowsInDiscretizedFramesR)
        
        flipIndexR = validRowsInDiscretizedFramesR.FlipIndex(iDiscFramesR);
        
        try
            % Find the corresponding Flip Index in flipDetailsL
            flipIndexL = find(processedPanelOutput.flipDetailsL.CorrectedTime == processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR));
        catch
        end

        % Check if a match was not found
        if(isempty(flipIndexL))
            % If no match is found, skip to the next iteration
            % disp(["Match wasn't found for " processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR) ", moving onto the next frame."])
            % val = [val;processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR)];
            continue;
        end
        
        % Find the index of the corresponding row in discretizedFramesL
        iDiscFramesL = find(discFramesL.FlipIndex == flipIndexL);
        
        % Assign the UnityMatchedFrame value from discretizedFramesR to discretizedFramesL
        discFramesL.UnityMatchedFrame(iDiscFramesL) = validRowsInDiscretizedFramesR.UnityMatchedFrame(iDiscFramesR);
    end
end

function [frameData] = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput)
    % Add a column with row indices to discretizedFramesR
    processedPanelOutput.discretizedFramesR.RowIndices = (1:height(processedPanelOutput.discretizedFramesR))';

    rowsWithAMatchedFrameR = ~isnan(processedPanelOutput.discretizedFramesR.UnityMatchedFrame);
    discretizedDataWithAMatchedFrameR = processedPanelOutput.discretizedFramesR(rowsWithAMatchedFrameR, :);

    rowsWithAMatchedFrameL = ~isnan(processedPanelOutput.discretizedFramesL.UnityMatchedFrame);
    discretizedDataWithAMatchedFrameL = processedPanelOutput.discretizedFramesL(rowsWithAMatchedFrameL, :);

    for iDiscRowR = 1:height(discretizedDataWithAMatchedFrameR)
        frameVal = discretizedDataWithAMatchedFrameR.UnityMatchedFrame(iDiscRowR);
        
        % Find the corresponding row in frameData
        matchingRow = frameData.Frame == frameVal;
        
        % Assign the original index to the frameData.DiscretizedFrameIndex
        frameData.DiscretizedFrameIndex(matchingRow) = discretizedDataWithAMatchedFrameR.RowIndices(iDiscRowR);

        % Assign FlipIndex to the frameData.DiscretizedFlipIndex 
        frameData.DiscretizedFlipIndex(matchingRow) = discretizedDataWithAMatchedFrameR.FlipIndex(iDiscRowR);
        
        % Assign the value of whether or not the flashed pattern matched
        % the expected pattern (Using the Right Panel)
        frameData.FrameValidity(matchingRow) = processedPanelOutput.frameDetailsR.Validity(iDiscRowR);
    end

    for iDiscRowL = 1:height(discretizedDataWithAMatchedFrameL)
        frameVal = discretizedDataWithAMatchedFrameL.UnityMatchedFrame(iDiscRowL);
        
        % Find the corresponding row in frameData
        matchingRow = frameData.Frame == frameVal;
        
        % Assign the Correct SynchBox Time to the
        % frameData.FrameOnsetSyncBoxTime (Using Left Panel)
        frameData.FrameOnsetSyncBoxTime(matchingRow) = processedPanelOutput.flipDetailsL.CorrectedTime(discretizedDataWithAMatchedFrameL.FlipIndex(iDiscRowL));
    end
end



