function frameSequenceIdxs = FindContinuousFrameSequences(frameDetails, frameVal)
    numRows = size(frameDetails, 1);
    frameSequenceIdxs = [];
    sequenceStart = 0; % Start index of a sequence of 0s
    
    for i = 1:numRows
        if frameDetails.Accuracy(i) == frameVal
            if sequenceStart == 0
                sequenceStart = i; % Start of a new sequence
            end
        else
            if sequenceStart > 0
                % End of the current sequence
                frameSequenceIdxs(end+1, :) = [sequenceStart, i - 1];
                sequenceStart = 0; % Reset for the next sequence
            end
        end
    end
    
    % Check if the last sequence goes until the end
    if sequenceStart > 0
        frameSequenceIdxs(end+1, :) = [sequenceStart, numRows];
    end
end

