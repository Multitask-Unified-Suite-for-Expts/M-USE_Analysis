function closestIndices = FindClosestWithWindow(TimeStamps, FlipTimes)
    % Initialize the output vector with zeros
    closestIndices = zeros(size(FlipTimes));
    % Length of the TimeStamps vector
    numTimeStamps = length(TimeStamps);
    
    % Define the window size
    windowSize = 1000;
    
    % Start index for the search window
    startIndex = 1;
    
    % Loop over each element in FlipTimes to find the closest TimeStamp
    for i = 1:length(FlipTimes)
        % Ending index for the search window
        endIndex = min(startIndex + windowSize - 1, numTimeStamps);
        
        % Get the current window from TimeStamps
        searchWindow = TimeStamps(startIndex:endIndex);
        
        % Calculate the absolute differences within the window
        differences = abs(searchWindow - FlipTimes(i));
        
        % Find the index of the minimum difference within the window
        [~, localIdx] = min(differences);
        
        % Calculate the global index
        closestIndices(i) = startIndex + localIdx - 1;
        
        % Update the start index for the next search window
        % Make sure the start index doesn't go beyond the length of TimeStamps
        startIndex = max(1, closestIndices(i) - windowSize / 2);
    end
end