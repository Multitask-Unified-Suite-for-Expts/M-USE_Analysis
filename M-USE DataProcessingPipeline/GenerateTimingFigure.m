function figureData = GenerateTimingFigure(serialRecvData)

    figureData = ProcessFlashPanelData(serialRecvData);

    DetectedFramesR = figureData.DetectedFramesR;

    unexpectedFrameSequenceIdxs = FindValidFrameSequences(DetectedFramesR, 0);
    expectedFrameSequenceIdxs = FindValidFrameSequences(DetectedFramesR, 1);

    %find all unexpected sequences with a specified width and buffer of
    %good data around
    unexpectedSequenceWidth = 1;
    unexpectedSequenceBuffer = 50;

    % validSequenceIndices = unexpectedFrameSequenceIdxs(GetValidSequences(unexpectedFrameSequenceIdxs, unexpectedSequenceWidth, unexpectedSequenceBuffer), 1);
    validSequenceIndices = GetValidSequences(unexpectedFrameSequenceIdxs, unexpectedSequenceWidth, unexpectedSequenceBuffer);

    %choose one sequence to plot
    plotIndex = unexpectedFrameSequenceIdxs(validSequenceIndices(1));

    plotIs = (plotIndex- unexpectedSequenceBuffer : plotIndex + unexpectedSequenceBuffer)';
    figureDetectedFrames = DetectedFramesR(plotIs, :);
    figureDetectedFrames.FrameDetailRow = plotIs;

    % rightOrderedFrameDetails = HandleFrameDetails(serialRecvData, figureData);
    % 
    % 
    % 
    % 
    % figureFrameDetails = [rightOrderedFrameDetails{1,68}; rightOrderedFrameDetails{1,69}; rightOrderedFrameDetails{1,70}];

    [figureAnalogData, figureFlipDetails] = convert_frame_details_to_analog_data(figureDetectedFrames, figureData);

    % Plot the analog data
    hold on
    plot(figureAnalogData.SynchBoxTime, figureAnalogData.LightR / 255);
    plot([44.3848000000000, 44.3848000000000], ylim, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Temporal Imprecision');

    scatter(figureFlipDetails.SynchBoxTime, figureFlipDetails.LightR / 255);
    
% Draw vertical gray lines at CorrectedTimes in filteredFlipDetails
% correctTimes = figureFlipDetails.CorrectedTime;
% for i = 1:length(correctTimes)
%     plot([correctTimes(i), correctTimes(i)], ylim, 'Color', [0, 0, 0]);
% end

% Assuming figureDiscretizedData is a table with two columns named Col1 and Col2
[~, idx, idxUnique] = unique(figureDetectedFrames.FlipIndex(:), 'stable');
uniqueDiscretizedData = figureDetectedFrames(idx, :);

% Create rectangles based on uniqueDiscretizedData
for i = 1:height(uniqueDiscretizedData)
    lightValue = uniqueDiscretizedData.Status(i);
    correspondingFlipRow = uniqueDiscretizedData.FlipIndex(i);
    correspondingAnalogIndex = figureData.DetectedFlipsR.AnalogIndex(correspondingFlipRow);
    correspondingSynchTime = figureData.analog_orig.SynchBoxTime(correspondingAnalogIndex);

    if lightValue == 0 && i < height(uniqueDiscretizedData)

        correspondingAnalogIndex2 = figureData.DetectedFlipsR.AnalogIndex(correspondingFlipRow + 1);
        correspondingSynchTime2 = figureData.analog_orig.SynchBoxTime(correspondingAnalogIndex2);
        % Specify the width of the rectangle
        rectangleWidth = correspondingSynchTime2 - correspondingSynchTime;
        
        % Specify the height of the rectangle
        rectangleHeight = 1;
    
        % Create a rectangle at the y of lightValue and the x of the correspondingSynchTime
        rectangle('Position', [correspondingSynchTime lightValue, rectangleWidth, rectangleHeight], 'FaceColor', [0, 0, 0, 0.2], 'EdgeColor', 'none');
    end
end

hold off


end


function [filteredAnalogData, filteredFlipDetails] = convert_frame_details_to_analog_data(detectedFrames, figureData)
    % filteredDiscretizedData = frame_details_to_discretized_data(filteredFrameDetails, figureData.DetectedFramesR);
    filteredFlipDetails = discretized_data_to_flip_details(detectedFrames, figureData.DetectedFlipsR);
    filteredAnalogData = figureData.analog_orig(filteredFlipDetails.AnalogIndex,:);
    
end

function filteredFlipDetails = discretized_data_to_flip_details(discretizedData, flipDetails)
    startFlipDetailsIdx = discretizedData.FlipIndex(1);
    endFlipDetailsIdx = discretizedData.FlipIndex(height(discretizedData));
    
    filteredFlipDetails = flipDetails(startFlipDetailsIdx:endFlipDetailsIdx, :);
end

% function filteredDiscretizedData = frame_details_to_discretized_data(frameDetails, discretizedData)
%     startIndex = frameDetails{1, 4};
%     endIndex = frameDetails{height(frameDetails), 4};
%     filteredDiscretizedData = discretizedData(startIndex:endIndex, :);
% end
