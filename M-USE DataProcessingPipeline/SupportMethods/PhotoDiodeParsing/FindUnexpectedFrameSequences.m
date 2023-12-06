function unexpectedFrameSequenceIdxs = FindUnexpectedFrameSequences(frameDetails)
    numRows = size(frameDetails, 1);
    unexpectedFrameSequenceIdxs = [];
    sequenceStart = 0; % Start index of a sequence of 0s
    
    for i = 1:numRows
        if frameDetails.Accuracy(i) == 0
            if sequenceStart == 0
                sequenceStart = i; % Start of a new sequence
            end
        else
            if sequenceStart > 0
                % End of the current sequence
                unexpectedFrameSequenceIdxs(end+1, :) = [sequenceStart, i - 1];
                sequenceStart = 0; % Reset for the next sequence
            end
        end
    end
    
    % Check if the last sequence goes until the end
    if sequenceStart > 0
        unexpectedFrameSequenceIdxs(end+1, :) = [sequenceStart, numRows];
    end
end

