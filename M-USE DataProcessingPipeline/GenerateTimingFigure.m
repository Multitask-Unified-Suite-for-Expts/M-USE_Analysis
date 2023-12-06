function figureData = GenerateTimingFigure(serialRecvData)

    figureData = ProcessFlashPanelData(serialRecvData);

    rightOrderedFrameDetails = HandleFrameDetails(serialRecvData, figureData.frameDetailsR);
    figureFrameDetails = [rightOrderedFrameDetails{1,68}; rightOrderedFrameDetails{1,69}; rightOrderedFrameDetails{1,70}];

    [figureAnalogData, figureFlipDetails, figureDiscretizedData] = convert_frame_details_to_analog_data(figureFrameDetails, figureData);

    % Assuming figureFrameDetails is a table with a column named Accuracy
    % Assuming figureFrameDetails is a table with a column named Accuracy
badFrameIdx = find(figureFrameDetails(:, 2) == 0);
badFrameDiscretizedIdx = figureFrameDetails(badFrameIdx, 4);
badFrameFlipIdx = figureData.discretizedFramesR(badFrameDiscretizedIdx, 2);
badFrameAnalogIdx = figureData.flipDetailsR(badFrameFlipIdx, :);

badFrameAnalogData = figureData.analog_orig(badFrameAnalogIdx, :);
    
    % Plot the analog data
    hold on
    plot(figureAnalogData.SynchBoxTime, figureAnalogData.LightR / 255);
    scatter(figureFlipDetails.SynchBoxTime, figureFlipDetails.LightR / 255);
    
% Draw vertical gray lines at CorrectedTimes in filteredFlipDetails
correctTimes = figureFlipDetails.CorrectedTime;
for i = 1:length(correctTimes)
    plot([correctTimes(i), correctTimes(i)], ylim, 'Color', [0, 0, 0]);
end

% Assuming figureDiscretizedData is a table with two columns named Col1 and Col2
[~, idx, idxUnique] = unique(figureDiscretizedData(:, 2), 'stable');
uniqueDiscretizedData = figureDiscretizedData(idx, :);

% Create rectangles based on uniqueDiscretizedData
for i = 1:height(uniqueDiscretizedData)
    lightValue = uniqueDiscretizedData{i, 1};
    correspondingFlipRow = uniqueDiscretizedData{i, 2};
    correspondingAnalogIndex = figureData.flipDetailsR.AnalogIndex(correspondingFlipRow);
    correspondingSynchTime = figureData.analog_orig.SynchBoxTime(correspondingAnalogIndex);

    % Specify the width of the rectangle
    rectangleWidth = 0.02;
    
    % Specify the height of the rectangle
    rectangleHeight = 0.05;

    % Create a rectangle at the y of lightValue and the x of the correspondingSynchTime
    rectangle('Position', [correspondingSynchTime lightValue - rectangleHeight/2, rectangleWidth, rectangleHeight], 'FaceColor', [0.8, 0.8, 0.8, 0.5], 'EdgeColor', 'none');
end

hold off


end


function [filteredAnalogData, filteredFlipDetails, filteredDiscretizedData] = convert_frame_details_to_analog_data(filteredFrameDetails, figureData)
    filteredDiscretizedData = frame_details_to_discretized_data(filteredFrameDetails, figureData.discretizedFramesR);
    filteredFlipDetails = discretized_data_to_flip_details(filteredDiscretizedData, figureData.flipDetailsR);
    filteredAnalogData = figureData.analog_orig(filteredFlipDetails.AnalogIndex,:);
    
end

function filteredDiscretizedData = frame_details_to_discretized_data(frameDetails, discretizedData)
    startIndex = frameDetails(1, 4);
    endIndex = frameDetails(height(frameDetails), 4);
    filteredDiscretizedData = discretizedData(startIndex:endIndex, :);
end

function filteredFlipDetails = discretized_data_to_flip_details(discretizedData, flipDetails)
    startFlipDetailsIdx = discretizedData{1,2};
    endFlipDetailsIdx = discretizedData{height(discretizedData),2};
    
    filteredFlipDetails = flipDetails(startFlipDetailsIdx:endFlipDetailsIdx, :);
end
