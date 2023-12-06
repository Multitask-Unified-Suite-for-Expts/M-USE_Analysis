function [accurateReportings] = HandleFrameDetails(serialRecvData, figureData)
    % figureData = ProcessFlashPanelData(serialRecvData);
    
    frameDetails = figureData.frameDetailsR;
    % [~, idx, ~] = unique(frameData.Frame, 'stable');
    % duplicates = frameData(ismember(1:height(frameData), idx), :);
    % 
    % newFrameData = duplicates;
    %accurateReportings = ProcessFrameDetails_OnlyAccurate(frameDetails);
badFrameDetails = GetBadFrameDetails(frameDetails);

    [accurateReportings] = ProcessFrameDetails_Both(frameDetails);

    % Display the ordered matrices
    nonEmptyMatrices = ~cellfun('isempty', accurateReportings);
    for i = find(nonEmptyMatrices)
        disp(['Matrix ' num2str(i)]);
        disp(accurateReportings{i});
    end
end

function badFrameDetails = GetBadFrameDetails(frameDetails)
    % 
    % ibadFrameDetails = find(frameDetails(:,2) == 0);
    % badFrameDetails = [frameDetails(ibadFrameDetails) ibadFrameDetails];


    % Find indices where the second column transitions from 1 to 0 or 0 to 1
    transitionIndices = find(diff(frameDetails{:, 2}) ~= 0);
    firstBadFrame = find(frameDetails{:,2} == 0,1);
    badSegmentStarts = [firstBadFrame
    fred = 2;

end


function [orderedMatrices] = ProcessFrameDetails_OnlyAccurate(frameDetails)
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
    
    % % Find the first occurrence of each pattern in the FlashPanelRStatus column
    % patterns = cellfun(@(mat) mat(:, 1), orderedMatrices, 'UniformOutput', false);
    % 
    % % Initialize array to store frame counts for each row in the orderedMatrices
    % frameCounts = zeros(size(appendedValues, 1), 1);
    % startIndex = 1;
    % % Search for each pattern in frameData.FlashPanelRStatus
    % for i = 1:numel(patterns)
    %     pattern = patterns{i};
    %     patternLength = numel(pattern);
    % 
    %     % Loop through frameData.FlashPanelRStatus to find the pattern
    %     for j = startIndex:(numel(frameData.FlashPanelRStatus) - patternLength + 1)
    %         if all(frameData.FlashPanelRStatus(j:(j + patternLength - 1)) == pattern)
    %             % Record the frame count for each row in the orderedMatrices
    %             frameCounts(i) = frameData.Frame(j);
    %             startIndex = j;
    %             break;
    %         end
    %     end
    % end
    
end

function [orderedMatrices] = ProcessFrameDetails_Both(frameDetails)
    % Ensure frameDetails is a table
    if ~istable(frameDetails)
        error('Input must be a table.');
    end
    
    % Find indices where the second column transitions from 1 to 0 or 0 to 1
    transitionIndices = find(diff(frameDetails{:, 2}) ~= 0);

    % Initialize cell array to store ordered matrices
    orderedMatrices = cell(1, numel(transitionIndices) + 1);
    
    % Initialize array to store appended values
    appendedValues = [];
    
    % Split the matrix into ordered matrices and include row numbers
    startIndex = 1;
    for i = 1:numel(transitionIndices)
        endIndex = transitionIndices(i);
        
        % Include all portions regardless of the second column value
        % Include the third column in the output matrix
        orderedMatrices{i} = [frameDetails{startIndex:endIndex, :}, (startIndex:endIndex)'];
        
        % Append new values to the array when switching from 0 to 1
        if endIndex + 1 <= size(frameDetails, 1) && frameDetails{endIndex + 1, 2} == 1
            % Include the third column in the appended values
            appendedValues = [appendedValues; frameDetails{endIndex + 1, :}, endIndex + 1];
        end
        
        startIndex = endIndex + 1;
    end
    
    % Include the last portion
    orderedMatrices{end} = [frameDetails{startIndex:end, :}, (startIndex:size(frameDetails, 1))'];
end
% 
%     % Find the first occurrence of each pattern in the FlashPanelRStatus column
%     patterns = cellfun(@(mat) mat(:, 1), orderedMatrices, 'UniformOutput', false);
% 
%     % Initialize array to store frame counts for each row in the orderedMatrices
%     frameCounts = zeros(size(appendedValues, 1), 1);
%     startIndex = 1;
% 
%     % Search for each pattern in frameData.FlashPanelRStatus
%     for i = 1:numel(patterns)
%         pattern = patterns{i};
%         patternLength = numel(pattern);
% 
%         % Loop through frameData.FlashPanelRStatus to find the pattern
%         for j = startIndex:(numel(frameData.FlashPanelRStatus) - patternLength + 1)
%             if all(frameData.FlashPanelRStatus(j:(j + patternLength - 1)) == pattern)
%                 % Record the frame count for each row in the orderedMatrices
%                 frameCounts(i) = frameData.Frame(j);
%                 startIndex = j;
%                 break;
%             end
%         end
%     end
% end

