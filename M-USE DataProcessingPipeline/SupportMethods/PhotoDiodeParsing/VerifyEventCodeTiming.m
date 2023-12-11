function [frameData, timeDifferences] = VerifyEventCodeTiming(frameData, serialRecvData, eventCode)
   
    % Find all the instances of the event code being sent in the frameData
    iFramesWithGivenEventCode = contains(frameData.FrameEventCodes, num2str(eventCode));
    
    % Use logical indexing to directly access frames with the specified event code
    framesWithGivenEventCode = frameData(iFramesWithGivenEventCode, :);
    timeDifferences = [];

    for idx = 1:height(framesWithGivenEventCode)
        frame = framesWithGivenEventCode.Frame(idx);
        refCode = str2double(framesWithGivenEventCode.ReferenceEventCodes{idx});
        detectedFrameOnset = framesWithGivenEventCode.DetectedFrameOnset(idx);
        % Find the first instance of the ref code in the
        % serialRecvData.EventCodes.CodeValue column where the
        % serialRecvData.EventCodes.CodeValue is >= frame
        if isnan(detectedFrameOnset)
            continue;
        end
        matchingIndex = find(serialRecvData.EventCodes.UnityRecvFrame >= frame & serialRecvData.EventCodes.CodeValue == refCode, 1);
        
        if ~isempty(matchingIndex)  
            % Process the found instance
            serialSynchBoxTime = serialRecvData.EventCodes.SynchBoxTime(matchingIndex);
            timeDifference = detectedFrameOnset - serialSynchBoxTime;
            timeDifferences = [timeDifferences; timeDifference];
            
            % Break the loop after processing the first match
            continue;
        else
            % Handle the case when no matching instance is found
            disp(['No matching instance found for frame ' num2str(frame) ' and reference code ' num2str(refCode)]);
        end
    end

    % Plot histogram
    figure;
    edges = -1:0.05:0.6;
    histogram(timeDifferences, edges, 'EdgeColor', 'black');
    title('Distribution of Time Differences');
    xlabel('Time Differences');
    ylabel('Frequency');
end
