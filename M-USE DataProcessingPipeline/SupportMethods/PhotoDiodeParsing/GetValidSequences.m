function validSequenceIndices = GetValidSequences(unexpectedFrameSequenceIdxs, sequenceWidth, sequenceBuffer)
    numSequences = size(unexpectedFrameSequenceIdxs, 1);
    validSequenceIndices = [];
    
    for i = 1:numSequences
        width = unexpectedFrameSequenceIdxs(i, 2) - unexpectedFrameSequenceIdxs(i, 1) + 1;
        
        if width == sequenceWidth
            if i == 1
                prevBuffer = true; % No previous sequence, assume sufficient buffer
            else
                prevBuffer = (unexpectedFrameSequenceIdxs(i, 1) - unexpectedFrameSequenceIdxs(i-1, 2)) >= sequenceBuffer;
            end
            
            if i == numSequences
                nextBuffer = true; % No next sequence, assume sufficient buffer
            else
                nextBuffer = (unexpectedFrameSequenceIdxs(i+1, 1) - unexpectedFrameSequenceIdxs(i, 2)) >= sequenceBuffer;
            end
            
            if prevBuffer && nextBuffer
                validSequenceIndices(end+1) = i;
            end
        end
    end
end