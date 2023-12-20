function [matchedFrameData, processedPanelOutput] = MatchSerialRecvAndFrameData(serialRecvData, frameData)
    processedPanelOutput = ProcessFlashPanelData(serialRecvData);
    processedPanelOutput.DetectedFramesR.UnityRecvFrame = processedPanelOutput.DetectedFlipsR.UnityRecvFrame(processedPanelOutput.DetectedFramesR.FlipIndex);
    processedPanelOutput.DetectedFramesL.UnityRecvFrame = processedPanelOutput.DetectedFlipsL.UnityRecvFrame(processedPanelOutput.DetectedFramesL.FlipIndex);

    
    processedPanelOutput.DetectedFramesR.ReportedMatchedFrame = nan(height(processedPanelOutput.DetectedFramesR),1);
    processedPanelOutput.DetectedFramesL.ReportedMatchedFrame = nan(height(processedPanelOutput.DetectedFramesL),1);
    % processedPanelOutput.DetectedFlipsR.ReportedMatchedOnsetFrame = nan(height(processedPanelOutput.DetectedFlipsR),1);
    % processedPanelOutput.DetectedFlipsL.ReportedMatchedOnsetFrame = nan(height(processedPanelOutput.DetectedFlipsL),1);
    % processedPanelOutput.DetectedFlipsR.DetectedMatchedOnsetFrame = nan(height(processedPanelOutput.DetectedFlipsR),1);
    % processedPanelOutput.DetectedFlipsL.DetectedMatchedOnsetFrame = nan(height(processedPanelOutput.DetectedFlipsL),1);
    frameData.FrameValidity = zeros(height(frameData), 1);    
    frameData.FrameOnsetSyncBoxTime = nan(height(frameData), 1);
    frameData.DetectedFrameIndex = nan(height(frameData), 1);
    frameData.FlipIndex = nan(height(frameData), 1);


    
    %Remove the rows that are showing duplicates
    frameData(diff(frameData.Frame) == 0, :) = [];
    % [~, idx, ~] = unique(frameData.Frame, 'stable');
    % nonduplicates = frameData(ismember(1:height(frameData), idx), :);
    % frameData = nonduplicates;

    validFrameSequenceIdxsR = FindValidFrameSequences(processedPanelOutput.DetectedFramesR, 1);
    for iSeq = 1:size(validFrameSequenceIdxsR,1)
        if validFrameSequenceIdxsR(iSeq,2) - validFrameSequenceIdxsR(iSeq,1) >= 24
            [processedPanelOutput.DetectedFramesR, frameData] = FindMatchingFrames(processedPanelOutput.DetectedFramesR, processedPanelOutput.DetectedFlipsR, validFrameSequenceIdxsR(iSeq,:), frameData, 'FlashPanelRStatus');
        else
            disp('asdkghfkdgh')
        end
    end

    [processedPanelOutput.DetectedFramesL] = AssignMatchedFramesToDetectedFramesLeft(processedPanelOutput);
    matchedFrameData = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput);
    fred = 2;
end

function [detectedFrames, unityReportedFrameData] = FindMatchingFrames(detectedFrames, detectedFlips, detectedFrameBoundaries, unityReportedFrameData, frameStatusCol)


    %find the frames in unityReportedFrameData that corresponed to the sequence in
    %detectedFrames between indices specified by detectedFrameBoundaries

    pad = 24;

    %exact window of detected frames
    detectedFrameSubset = detectedFrames(detectedFrameBoundaries(1):detectedFrameBoundaries(2),:);


    %flipDetails contains frames in which unity received syncbox data
    flipSubsetIdxs = max(detectedFrameSubset.FlipIndex(1) - pad, 1) : min(detectedFrameSubset.FlipIndex(end) + pad, height(detectedFlips));
    try
        detectedFlipSubset = detectedFlips(flipSubsetIdxs,:);
    catch
        fred = 2;
    end
    reportedFrameSubsetIdxs = find(unityReportedFrameData.Frame >= detectedFlipSubset.UnityRecvFrame(1) & unityReportedFrameData.Frame < detectedFlipSubset.UnityRecvFrame(end));
    reportedFrameSubset = unityReportedFrameData(reportedFrameSubsetIdxs,:);

    reportedFrameStatus = reportedFrameSubset.(frameStatusCol);
    detectedFrameStatus = detectedFrameSubset.Status;

    %find matches to the detected frame vector within the reported frames
    matchPoints = strfind(reportedFrameStatus', detectedFrameStatus'); 

    %the light sensor data from the first frame in the detected frame 
    % sequence was sent back to unity on a given frame, and the match point
    % cannot be greater than this
    maxMatchPoint = find(reportedFrameSubset.Frame == detectedFlips.UnityRecvFrame(detectedFrameSubset.FlipIndex(1)), 1, "last");

    matchPoints(matchPoints > maxMatchPoint) = [];

    if ~isempty(matchPoints)
        matchPoint = max(matchPoints); %get the closest point to the received frame
        reportedMatchIdxs = reportedFrameSubsetIdxs(1) + matchPoint - 1 : reportedFrameSubsetIdxs(1) + matchPoint - 1 + height(detectedFrameSubset) - 1;
        unityReportedFrameData.FrameValidity(reportedMatchIdxs) = 1;
        unityReportedFrameData.DetectedFrameIndex(reportedMatchIdxs) = detectedFrameBoundaries(1):detectedFrameBoundaries(2);
        unityReportedFrameData.FlipIndex(reportedMatchIdxs) = detectedFrameSubset.FlipIndex;
        detectedFrames.ReportedMatchedFrame(detectedFrameBoundaries(1):detectedFrameBoundaries(2)) = unityReportedFrameData.Frame(reportedMatchIdxs);
        detectedFlipMatchIdxs = flipSubsetIdxs(1) + matchPoint - 1 : flipSubsetIdxs(1) + matchPoint - 1 + height(detectedFrameSubset) - 1;
    else
        disp("OH NO!!!!");
    end

    

    % lag = 0;
    % reportedFrameStatus = originalReportedFrameStatus;
    % 
    % 
    % % Create a comparison matrix
    % try
    %     comparison = [reportedFrameStatus(1:length(detectedFrameStatus)), detectedFrameStatus];
    % catch
    %     fred = 2;
    %     comparison = [reportedFrameStatus, detectedFrameStatus(1:length(reportedFrameStatus))];
    % 
    % end
    % % Check for equality
    % if ~isequal(comparison(:, 1), comparison(:, 2))
    % 
    %     shiftFactor = 0;
    %     % sequenceIs(2) = sequenceIs(2) - 1;
    %     breakOccurred = false;  % Initialize the variable
    % 
    %     % Keep shifting until the sequences match or lag reaches a limit
    %     while ~isequal(comparison(:, 1), comparison(:, 2))
    %         %disp("SHIFTING")
    %         lag = lag + 1;
    % 
    %          if (lag > 24)
    %              breakOccurred = true;
    %              disp(["LAG " lag]);
    %              lag = 1;
    % 
    %              break;
    %          end
    % 
    %         % discFrameVals = DetectedFramesR(sequenceIs(1):(sequenceIs(2)-1), :);
    %         % flipDetailVals = flipDetails(max(discFrameVals{1, 2}, 1) : min(discFrameVals{end, 2} + 24, end), :);
    %         % 
    %         % frameVals = frameData(frameData.Frame >= flipDetailVals.UnityRecvFrame(1) & frameData.Frame <= flipDetailVals.UnityRecvFrame(end), :);
    % 
    %         % Update reportedFrameStatus after each shift
    %         try
    %             reportedFrameStatus = [originalReportedFrameStatus(lag:end); originalReportedFrameStatus(1:lag-1)];
    %         catch
    %             disp("out of range")
    %         end
    %         % Update comparison matrix
    %         comparison = [reportedFrameStatus(1:length(detectedFrameStatus)), detectedFrameStatus];
    %     end
    % 
    %     try
    %         detectedFrames.UnityMatchedFrame(detectedFrameBoundaries(1):(detectedFrameBoundaries(2)-1)) = reportedFrameSubset.Frame(lag:lag+length(detectedFrameStatus) - 1);
    %     catch
    %         detectedFrames.UnityMatchedFrame(detectedFrameBoundaries(1):detectedFrameBoundaries(1)+height(reportedFrameSubset.Frame(lag:end))-1) = reportedFrameSubset.Frame(lag:end);
    % 
    %         fred = 2;
    %     end
    % else
    %    % disp("EQUALITY")
    %    % Determine the valid range based on the lengths of sequence and frameVals.Frame
    % range = detectedFrameBoundaries(1):(detectedFrameBoundaries(1) + length(reportedFrameSubset.Frame(1:length(detectedFrameStatus))) - 1);
    % 
    % % Update only the existing rows in DetectedFramesR.UnityMatchedFrame
    % detectedFrames.UnityMatchedFrame(range) = reportedFrameSubset.Frame(1:length(detectedFrameStatus));
    % 
    % end
end

function [discFramesL] = AssignMatchedFramesToDetectedFramesLeft(processedPanelOutput)
    
    % Find the rows in the DetectedFramesR that have a matched Frame 
    validRowsInDetectedFramesR = processedPanelOutput.DetectedFramesR(~isnan(processedPanelOutput.DetectedFramesR.ReportedMatchedFrame), :);
    
    % Get a reference to the DetectedFramesL
    discFramesL = processedPanelOutput.DetectedFramesL;
    
    % Loop through each valid row in DetectedFramesR
    for iDiscFramesR = 1:height(validRowsInDetectedFramesR)
        
        flipIndexR = validRowsInDetectedFramesR.FlipIndex(iDiscFramesR);
        
        try
            % Find the corresponding Flip Index in flipDetailsL
            flipIndexL = find(processedPanelOutput.DetectedFlipsL.CorrectedTime == processedPanelOutput.DetectedFlipsR.CorrectedTime(flipIndexR));
        catch
        end

        % Check if a match was not found
        if(isempty(flipIndexL))
            % If no match is found, skip to the next iteration
            % disp(["Match wasn't found for " processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR) ", moving onto the next frame."])
            % val = [val;processedPanelOutput.flipDetailsR.CorrectedTime(flipIndexR)];
            continue;
        end
        
        % Find the index of the corresponding row in DetectedFramesL
        iDiscFramesL = find(discFramesL.FlipIndex == flipIndexL);
        
        % Assign the UnityMatchedFrame value from DetectedFramesR to DetectedFramesL
        discFramesL.ReportedMatchedFrame(iDiscFramesL) = validRowsInDetectedFramesR.ReportedMatchedFrame(iDiscFramesR);
    end
end

function [frameData] = AssignDiscretizedDataFieldsToFrameData(frameData, processedPanelOutput)
    % % Add a column with row indices to DetectedFramesR
    % processedPanelOutput.DetectedFramesR.RowIndices = (1:height(processedPanelOutput.DetectedFramesR))';
    % 
    % rowsWithAMatchedFrameR = ~isnan(processedPanelOutput.DetectedFramesR.UnityMatchedFrame);
    % discretizedDataWithAMatchedFrameR = processedPanelOutput.DetectedFramesR(rowsWithAMatchedFrameR, :);
    % 
    rowsWithAMatchedFrameL = ~isnan(processedPanelOutput.DetectedFramesL.ReportedMatchedFrame);
    discretizedDataWithAMatchedFrameL = processedPanelOutput.DetectedFramesL(rowsWithAMatchedFrameL, :);
    % 
    % for iDiscRowR = 1:height(discretizedDataWithAMatchedFrameR)
    %     frameVal = discretizedDataWithAMatchedFrameR.UnityMatchedFrame(iDiscRowR);
    % 
    %     % Find the corresponding row in frameData
    %     matchingRow = frameData.Frame == frameVal;
    % 
    %     % Assign the original index to the frameData.DiscretizedFrameIndex
    %     frameData.DiscretizedFrameIndex(matchingRow) = discretizedDataWithAMatchedFrameR.RowIndices(iDiscRowR);
    % 
    %     % Assign FlipIndex to the frameData.DiscretizedFlipIndex 
    %     frameData.DiscretizedFlipIndex(matchingRow) = discretizedDataWithAMatchedFrameR.FlipIndex(iDiscRowR);
    % 
    %     % Assign the value of whether or not the flashed pattern matched
    %     % the expected pattern (Using the Right Panel)
    %     frameData.FrameValidity(matchingRow) = processedPanelOutput.frameDetailsR.Validity(iDiscRowR);
    % end

    for iDiscRowL = 1:height(discretizedDataWithAMatchedFrameL)
        frameVal = discretizedDataWithAMatchedFrameL.ReportedMatchedFrame(iDiscRowL);
        
        % Find the corresponding row in frameData
        matchingRow = frameData.Frame == frameVal;
        
        % Assign the Correct SynchBox Time to the
        % frameData.FrameOnsetSyncBoxTime (Using Left Panel)
        frameData.FrameOnsetSyncBoxTime(matchingRow) = processedPanelOutput.DetectedFlipsL.CorrectedTime(discretizedDataWithAMatchedFrameL.FlipIndex(iDiscRowL));
    end
end



