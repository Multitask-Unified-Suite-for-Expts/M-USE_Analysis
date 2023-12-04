function problemFrames = DetectPatternViolations(discretizedFrames, expectedPattern, comparisonWindowSize)


% Preallocate the matrix
allPossibleExpectedPatterns = zeros(length(expectedPattern), length(expectedPattern));
% Loop to fill in the matrix
for i = 1:length(expectedPattern)
    allPossibleExpectedPatterns(:, i) = circshift(expectedPattern, i - 1);
end

numPatterns = ceil(comparisonWindowSize / length(expectedPattern));
patternedWindow = repmat(expectedPattern, 1, numPatterns);

numWindows = ceil(length(discretizedFrames) / length(patternedWindow));

problemFrames = zeros(size(discretizedFrames));


%find first exact match to expectedPattern
lag = FindFirstMatch(discretizedFrames, expectedPattern);
if lag > 1
    problemFrames = (1:lag-1)';
end


problemWindows = [];
for iWin = 1 : numWindows

    detectedVals = discretizedFrames(iWin * comparisonWindowSize - comparisonWindowSize + 1 : min(iWin * comparisonWindowSize, end));

    r = corr(detectedVals, patternedWindow(1:length(detectedVals))');
    if r < 1
        firstMatch = FindFirstMatch(detectedVals, expectedPattern);

        problemWindows = [problemWindows; iWin];
        firstViolation = find(detectedVals ~= patternedWindow(1:length(detectedVals))', 1);
        sum(abs(detectedVals - patternedWindow(1:length(detectedVals))'));
        disp('alskhglkasdhg');
    end
    %need to recursively search for exceptions, keep halving window size
    %when results are not perfect

end
problemFrames = 1;

end

function lag = FindFirstMatch(detectedData, expectedPattern)
%find first exact match to expectedPattern
match = 0;
lag = 0;
while match == 0
    comparisonData = discretizedFrames(lag + 1 : length(expectedPattern) + lag);

    if(size(expectedPattern, 1) == 1 && size(comparisonData, 1) ~= 1)
        comparisonData = comparisonData';
    end

    if isequal(comparisonData, expectedPattern)
        match = 1;
    end
    
    lag = lag + 1;
end
end



