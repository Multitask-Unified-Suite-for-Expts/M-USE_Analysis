function [frameDetails, exceptionDetails] = FlashPanelPatternProcessor(discretizedFrames, expectedPattern)

    % Matrix of all possible orders of the pattern
    patternTemplate = zeros(length(expectedPattern), length(expectedPattern));
    % Loop to fill in the matrix
    for i = 1:length(expectedPattern)
        patternTemplate(:, i) = circshift(expectedPattern, i - 1);
    end
    
    patternMatchers = GeneratePatternMatchers(patternTemplate, length(expectedPattern));
    
    iFrame = 1;
    patternCol = 1;
    exceptionFound = 1;

    frameDetails = [discretizedFrames zeros(length(discretizedFrames),1) zeros(length(discretizedFrames),1)];
    exceptionDetails.NoMatch = [];
    exceptionDetails.Stutter = [];
    exceptionDetails.PatternBreak = [];
    
    while iFrame < length(discretizedFrames)
        if exceptionFound
            [iFrame, patternCol, frameDetails, exceptionDetails] = FindNextMatch(discretizedFrames, patternTemplate, iFrame, patternCol, frameDetails, exceptionDetails);
        end
        if iFrame < length(discretizedFrames)
            [iFrame, frameDetails, exceptionDetails, exceptionFound] = FindNextException(discretizedFrames, patternMatchers, iFrame, patternCol, frameDetails, exceptionDetails);
        end
    end

end


function patternMatchers = GeneratePatternMatchers(patternTemplate, maxSize)

    numPatterns = floor(maxSize / size(patternTemplate,1));
    patternMatchers{1} = repmat(patternTemplate, 1, numPatterns);
    iP = 2;

    while numPatterns >= 2
        numPatterns = floor(numPatterns / 2);
        patternMatchers{iP} = repmat(patternTemplate, 1, numPatterns);
        iP = iP + 1;
    end

end

function [iFrame, patternCol, frameDetails, exceptionDetails] = FindNextMatch(discretizedFrames, patternTemplate, iFrame, patternCol, frameDetails, exceptionDetails)

    dFrames = discretizedFrames(iFrame:end);
    iFr = 1;
    matchCol = 0;

    while matchCol == 0

        %find if current segment matches any of the columns in
        %patternTemplate
        if iFr <= (size(dFrames,1) - size(patternTemplate,1)) + 1
            matchCol = FindMatchColumn(patternTemplate, dFrames, iFr);
        elseif iFr <= size(dFrames,1)
            ht = size(dFrames,1) - iFr + 1;
            %there might be more than one match to patterns shorter than
            %the template, so loop around from the currently selected pattern
            %column to
            matchCol = FindMatchColumn(patternTemplate(1:ht,patternCol:end), dFrames, iFr, 1);
            if isempty(matchCol)
                matchCol = FindMatchColumn(patternTemplate(1:ht,1:patternCol-1), dFrames, iFr);
            end
        else
            matchCol = -1;
        end

        if isempty(matchCol) 
            %no matches, move to next frame
            matchCol = 0;
            exceptionDetails.NoMatch = [exceptionDetails.NoMatch; iFrame + iFr - 1 0];
            frameDetails(iFrame + iFr, 2:3) = [0 matchCol];
            iFr = iFr + 1;
        else
            %we found the match
            iFrame = iFrame + iFr - 1; 
            frameDetails(iFrame, 2:3) = [1 matchCol];
        end
    end

    if matchCol ~= patternCol
        exceptionDetails.Stutter = [exceptionDetails.Stutter; iFrame matchCol];
            frameDetails(iFrame, 2:3) = [0 matchCol];
    end

    patternCol = matchCol;

end

function matchCol = FindMatchColumn(patternTemplate, dFrames, iFr, varargin)

    if length(varargin) == 1 && varargin{1} == 1
        matchCol = find(sum(abs(patternTemplate - dFrames(iFr:iFr + size(patternTemplate,1) - 1))) == 0, 1);
    else
        matchCol = find(sum(abs(patternTemplate - dFrames(iFr:iFr + size(patternTemplate,1) - 1))) == 0);
    end
    if length(matchCol) > 1
        disp('ERROR - more than one column match identified')
    end
    
end


function [iFrame, frameDetails, exceptionDetails, exceptionFound] = FindNextException(discretizedFrames, patternMatchers, iFrame, patternCol, frameDetails, exceptionDetails)

    dFrames = discretizedFrames(iFrame:end);
    
    exceptionFound = 1;
    for iP = 1:length(patternMatchers)
        pm = patternMatchers{iP}(:,patternCol);
        if length(pm) <= length(dFrames)
            if isequal(pm, dFrames(1:length(pm))) %perfect match
                exceptionFound = 0;
                frameDetails(iFrame + 1 : iFrame + length(pm), 2:3) = repmat([1 patternCol], length(pm),1);
                iFrame = iFrame + length(pm);
                break;
            end
            %if it's not a perfect match, go to the next half length pattern and
            %check
        end
    end

    if iP == length(patternMatchers) && exceptionFound %there is an exception to the smallest complete pattern 
        
        while length(pm) > 1
            pm = pm(1:end-1);
            if isequal(pm(1:min([length(pm) length(dFrames)])), dFrames(1:min([length(pm) length(dFrames)]))) %perfect match
                frameDetails(iFrame + 1 : iFrame + length(pm), 2:3) = repmat([1 patternCol], length(pm),1);
                iFrame = iFrame + length(pm);
                exceptionDetails.PatternBreak = [exceptionDetails.PatternBreak; iFrame patternCol];
                break;
            elseif length(pm) == 1

                frameDetails(iFrame + 1, 2:3) = [0 patternCol];
                iFrame = iFrame + 1;
                exceptionDetails.PatternBreak = [exceptionDetails.PatternBreak; iFrame patternCol];
            end
        end
    end

end