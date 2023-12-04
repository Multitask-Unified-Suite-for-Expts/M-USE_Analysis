function [frameDetailsL, exceptionDetailsL, frameDetailsR, exceptionDetailsR, accurateReportings, frameCounts] = HandleFrameDetails(serialRecvData, frameData)
    [frameDetailsL, exceptionDetailsL, frameDetailsR, exceptionDetailsR] = ProcessFlashPanelData(serialRecvData);
    
    % Process frameDetailsL
    frameDetails = frameDetailsR;
    [accurateReportings, frameCounts] = ProcessFrameDetails(frameDetails, frameData);

    % Display the ordered matrices
    nonEmptyMatrices = ~cellfun('isempty', accurateReportings);
    for i = find(nonEmptyMatrices)
        disp(['Matrix ' num2str(i)]);
        disp(accurateReportings{i});
    end
end

function [orderedMatrices, frameCounts] = ProcessFrameDetails(frameDetails, frameData)
    % Find indices where the second column transitions from 1 to 0
    splitIndices = find(diff([0; frameDetails(:, 2); 0]) == -1);
    
    % Initialize cell array to store ordered matrices
    orderedMatrices = cell(1, numel(splitIndices));
    
    % Initialize array to store appended values
    appendedValues = [];
    
    % Split the matrix into ordered matrices and include row numbers
    startIndex = 1;
    for i = 1:numel(splitIndices)
        endIndex = splitIndices(i) - 1;
        
        % Only include non-empty portions where the second column is 1
        if any(frameDetails(startIndex:endIndex, 2) == 1)
            orderedMatrices{i} = [frameDetails(startIndex:endIndex, :), (startIndex:endIndex)'];
        end
        
        % Append new values to the array when switching from 0 to 1
        if endIndex + 2 <= size(frameDetails, 1) && frameDetails(endIndex + 2, 2) == 1
            appendedValues = [appendedValues; frameDetails(endIndex + 2, :), endIndex + 2];
        end
        
        startIndex = endIndex + 2;
    end
    
    % Filter the matrices to include only rows where the second column is 1
    orderedMatrices = cellfun(@(mat) mat(mat(:, 2) == 1, :), orderedMatrices, 'UniformOutput', false);
    
    % Find the first occurrence of each pattern in the FlashPanelRStatus column
    patterns = cellfun(@(mat) mat(:, 1), orderedMatrices, 'UniformOutput', false);
    
    % Initialize array to store frame counts for each row in the orderedMatrices
    frameCounts = zeros(size(appendedValues, 1), 1);
    
    % Search for each pattern in frameData.FlashPanelRStatus
    for i = 1:numel(patterns)
        pattern = patterns{i};
        patternLength = numel(pattern);
        
        % Loop through frameData.FlashPanelRStatus to find the pattern
        for j = 1:(numel(frameData.FlashPanelRStatus) - patternLength + 1)
            if all(frameData.FlashPanelRStatus(j:(j + patternLength - 1)) == pattern)
                % Record the frame count for each row in the orderedMatrices
                frameCounts(i) = frameData.Frame(j);
                break;
            end
        end
    end
    
end
