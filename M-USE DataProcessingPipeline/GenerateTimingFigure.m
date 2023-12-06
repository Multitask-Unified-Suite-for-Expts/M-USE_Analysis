function figureData = GenerateTimingFigure(serialRecvData)

    figureData = ProcessFlashPanelData(serialRecvData);

    frameDetailsR = figureData.frameDetailsR;

    unexpectedFrameSequenceIdxs = FindContinuousFrameSequences(frameDetailsR, 0);
    expectedFrameSequenceIdxs = FindContinuousFrameSequences(frameDetailsR, 1);

    %find all unexpected sequences with a specified width and buffer of
    %good data around
    unexpectedSequenceWidth = 3;
    unexpectedSequenceBuffer = 50;

    validSequenceIndices = unexpectedFrameSequenceIdxs(GetValidSequences(unexpectedFrameSequenceIdxs, unexpectedSequenceWidth, unexpectedSequenceBuffer), 1);

    %choose one sequence to plot
    plotIndex = validSequenceIndices(1);

    plotIs = (plotIndex- unexpectedSequenceBuffer : plotIndex + unexpectedSequenceBuffer)';
    figureFrameDetails = frameDetailsR(plotIs, :);
    figureFrameDetails.FrameDetailRow = plotIs;

    % rightOrderedFrameDetails = HandleFrameDetails(serialRecvData, figureData);
    % 
    % 
    % 
    % 
    % figureFrameDetails = [rightOrderedFrameDetails{1,68}; rightOrderedFrameDetails{1,69}; rightOrderedFrameDetails{1,70}];

    [figureAnalogData, figureFlipDetails, figureDiscretizedData] = convert_frame_details_to_analog_data(figureFrameDetails, figureData);

    % Plot the analog data
    hold on
    plot(figureAnalogData.SynchBoxTime, figureAnalogData.LightR / 255);
    scatter(figureFlipDetails.SynchBoxTime, figureFlipDetails.LightR / 255);
    
% Draw vertical gray lines at CorrectedTimes in filteredFlipDetails
% correctTimes = figureFlipDetails.CorrectedTime;
% for i = 1:length(correctTimes)
%     plot([correctTimes(i), correctTimes(i)], ylim, 'Color', [0, 0, 0]);
% end

% Assuming figureDiscretizedData is a table with two columns named Col1 and Col2
[~, idx, idxUnique] = unique(figureDiscretizedData(:, 2), 'stable');
uniqueDiscretizedData = figureDiscretizedData(idx, :);

% Create rectangles based on uniqueDiscretizedData
for i = 1:height(uniqueDiscretizedData)
    lightValue = uniqueDiscretizedData{i, 1};
    correspondingFlipRow = uniqueDiscretizedData{i, 2};
    correspondingAnalogIndex = figureData.flipDetailsR.AnalogIndex(correspondingFlipRow);
    correspondingSynchTime = figureData.analog_orig.SynchBoxTime(correspondingAnalogIndex);

    if lightValue == 0 && i < height(uniqueDiscretizedData)

        correspondingAnalogIndex2 = figureData.flipDetailsR.AnalogIndex(correspondingFlipRow + 1);
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


function [filteredAnalogData, filteredFlipDetails, filteredDiscretizedData] = convert_frame_details_to_analog_data(filteredFrameDetails, figureData)
    filteredDiscretizedData = frame_details_to_discretized_data(filteredFrameDetails, figureData.discretizedFramesR);
    filteredFlipDetails = discretized_data_to_flip_details(filteredDiscretizedData, figureData.flipDetailsR);
    filteredAnalogData = figureData.analog_orig(filteredFlipDetails.AnalogIndex,:);
    
end

function filteredDiscretizedData = frame_details_to_discretized_data(frameDetails, discretizedData)
    startIndex = frameDetails{1, 4};
    endIndex = frameDetails{height(frameDetails), 4};
    filteredDiscretizedData = discretizedData(startIndex:endIndex, :);
end

function filteredFlipDetails = discretized_data_to_flip_details(discretizedData, flipDetails)
    startFlipDetailsIdx = discretizedData{1,2};
    endFlipDetailsIdx = discretizedData{height(discretizedData),2};
    
    filteredFlipDetails = flipDetails(startFlipDetailsIdx:endFlipDetailsIdx, :);
end
