function Import_MUSE_Session_RawData(varargin)
    % Set the current folder and configure script cleanup
    scriptPath = pwd;
    finishup = onCleanup(@() CleanupFun(scriptPath));
    % Optionally, add the current path to the MATLAB search path.
    % addpath(genpath(scriptPath));

    %% Set Subject Data Paths and Read Data
    % Define the experiment type or use 'FLU' as the default
    exptType = CheckVararginPairs('exptType', 'FLU', varargin{:});

    % Specify the data folder or prompt the user to choose one
    dataFolder = CheckVararginPairs('dataFolder', '', varargin{:});
    % folderName = CheckVararginPairs('folderName', '', varargin{:}); % Unused variable
    if isempty(dataFolder)
        dataFolder = uigetdir('Data', 'Choose the data folder');
    end

    % Specify gaze arguments
    gazeArgs = CheckVararginPairs('gazeArgs', '', varargin{:});

    %% Variables to Control Recalculation
    forceProcessAllData = 0;
    forceGazeClassification = 0;
    forceGazeFrameAlignment = 0;
    forceSaveFrameDataForReplayer = 0;
    ignoreGaze = 0;
    ignoreGazeClassification = 0;
    ignoreBlockCondition = 1;
    ignoreReplayer = 1;
    forceLP = 0;
    syncboxSerialData = 1;

    deleteAbortedTrials = 0;
    nanAbortedTrials = 1;

    % Find task names in the specified data folder
    taskNames = findTaskNames(dataFolder);
    taskSelectionNames = cellfun(@(x) x(1:8), taskNames, 'UniformOutput', false);
    taskSelectionNames = cellfun(@(name) ['TaskSelectionData/' name], taskSelectionNames, 'UniformOutput', false);

    % Define subject information
    subjectID = CheckVararginPairs('SubjectID', 'Subject1', varargin{:});
    subjectNum = CheckVararginPairs('Subject#', 1, varargin{:});

    sessionID = CheckVararginPairs('SessionID', 'Session1', varargin{:});
    sessionNum = CheckVararginPairs('Session#', 1, varargin{:});

    % Display a header for the processing step
    disp('-----------------------------------------------------------------')
    fprintf(['\n\nProcessing data files for subject ' subjectID ', session ' sessionID '.\n']);

    % Define runtime and processed data paths for each task
    [runTimeDataPaths, processedDataPaths] = defineDataStoragePaths(subjectID, dataFolder, taskNames);
    [taskSelectionRunTimeDataPaths, taskSelectionProcessedDataPaths] = defineDataStoragePaths(subjectID, dataFolder, taskSelectionNames);
    % Add "TaskSelection/" to each element using cellfun

    % Begin processing the data for each of the task folders
    for i = 1:numel(taskNames)
        currentTaskName = taskNames{i};
        currentTaskDataPath = runTimeDataPaths{i};
        currentProcessedDataPath = processedDataPaths{i};
        
        processAllDataTypes(currentTaskDataPath, currentTaskName, forceProcessAllData,  currentProcessedDataPath, ...
            subjectNum, sessionNum, deleteAbortedTrials, nanAbortedTrials, exptType, ...
            gazeArgs, ignoreBlockCondition);
        
        currentTaskSelectionFolderName = taskSelectionNames{i};
        currentTaskSelectionDataPath = taskSelectionRunTimeDataPaths{i};
        currentTaskSelectionProcessedDataPath = taskSelectionProcessedDataPaths{i};

         processTaskSelectionDataTypes(currentTaskSelectionDataPath, currentTaskSelectionFolderName, forceProcessAllData,  currentTaskSelectionProcessedDataPath, ...
            subjectNum, sessionNum, deleteAbortedTrials, nanAbortedTrials, exptType, ...
            gazeArgs, ignoreBlockCondition);
    end

    % Begin processing the data for each of the folders in task selection
    
end

function [runtimeDataPaths, processedDataPaths] = defineDataStoragePaths(subjectID, dataFolder, taskNames)
    global singleDataFileName singleDataFileVarNames multiDataFileVarNames; 
    
    % Initialize cell arrays to store runtime and processed data paths
    runtimeDataPaths = cell(numel(taskNames), 1);
    processedDataPaths = cell(numel(taskNames), 1);

    % Loop through each task to define data paths
    for i = 1:numel(taskNames)
        currentTaskName = taskNames{i};
        runtimeDataPaths{i} = [dataFolder filesep currentTaskName];
        processedDataPaths{i} = [dataFolder filesep 'ProcessedData' filesep currentTaskName];
        
        % Define the single data file name based on subject ID
        singleDataFileName = [subjectID '__ProcessedData.mat'];
        processedDataFilePath = [processedDataPaths{i} filesep singleDataFileName];

        % Create the processed data folder if it doesn't exist
        if ~exist(processedDataPaths{i}, 'dir')
            mkdir(processedDataPaths{i});
        end

        % Check if the single data file exists and get its variable names
        if exist(processedDataFilePath, 'file')
            fprintf('\tComplete subject data file exists, determining which variables are in it.\n');
            singleDataFileVarNames = who('-file', processedDataFilePath);
        else
            singleDataFileName = '';  % No single data file found
            singleDataFileVarNames = '';  % No variable names available
        end

        % Get the names of multi-data files in the processed data folder
        multiDataFileVarNames = dir([processedDataPaths{i} filesep '*.mat']);
        multiDataFileVarNames = {multiDataFileVarNames.name};
    end
end


function processAllDataTypes(runtimeDataPath, currentTaskName, forceProcessAllData, currentProcessedDataPath, ...
    subjectNum, sessionNum, deleteAbortedTrials, nanAbortedTrials, exptType, gazeArgs, ...
    ignoreBlockCondition)

    fprintf(['\n\nProcessing data for task ' currentTaskName '.\n']);
    
    % Load or process trial data
    trialData = processTrialData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, deleteAbortedTrials, ...
        nanAbortedTrials, subjectNum, sessionNum);

    % Load or process block data
    blockData = processBlockData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum, exptType, ignoreBlockCondition);

    % Load or process frame data
    frameData = processFrameData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum);

    % Load or process gaze data if gazeArgs is not 'cancel'
    if ~strcmp(gazeArgs, 'cancel')
        gazeData = processGazeData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum);
    else
        gazeData = [];
    end
    
    % Load or process serial data
    processSerialData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum);
end

function processTaskSelectionDataTypes(runtimeDataPath, currentTaskSelectionFolderName, forceProcessAllData, currentProcessedDataPath, ...
    subjectNum, sessionNum, deleteAbortedTrials, nanAbortedTrials, exptType, gazeArgs, ...
    ignoreBlockCondition)

    fprintf(['\n\nProcessing data for task ' currentTaskSelectionFolderName '.\n']);

  % Load or process frame data
    frameData = processFrameData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum);

    % Load or process gaze data if gazeArgs is not 'cancel'
    if ~strcmp(gazeArgs, 'cancel')
        gazeData = processGazeData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum);
    else
        gazeData = [];
    end
    
    % Load or process serial data
    [serialSentData, serialRecvData] = processSerialData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum);
end


%% Process Trial Data
function trialData = processTrialData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, deleteAbortedTrials, nanAbortedTrials, subjectNum, sessionNum)
    % Display a message indicating that trial data is being loaded.
    fprintf('\tLoading trial data from file.\n');
    
    % Attempt to load trial data from a processed data file.
    trialData = LoadDataCheckCompleteFile('trialData', currentProcessedDataPath, forceProcessAllData);

    % If trial data is not loaded, read it from raw data files.
    if isempty(trialData)
        % Read trial data from raw data files.
        trialData = ReadDataFiles([runtimeDataPath filesep 'TrialData'], '*TrialData.txt', 'ImportOptions', {'delimiter', '\t', 'TreatAsEmpty',{'null'}});
        
        % Find aborted trials based on the 'AbortCode' column.
        abortedTrials = find(trialData.AbortCode > 0);
        abortedTrialData = trialData(abortedTrials, :);

        % Handle aborted trials based on configuration options.
        if deleteAbortedTrials
            trialData(abortedTrials, :) = [];
        elseif nanAbortedTrials
            vars = trialData.Properties.VariableNames;
            abortCol = find(strcmp(vars, 'AbortCode'));

            % Replace aborted trial values with empty or NaN.
            for i = abortedTrials'
                for iCol = abortCol + 1 : width(trialData)
                    if iscell(trialData.(vars{iCol}))
                        trialData.(vars{iCol}){i} = '';
                    elseif isnumeric(trialData.(vars{iCol}))
                        trialData.(vars{iCol})(i) = NaN;
                    end
                end
            end
        end

        % Add subject and session information to trial data.
        trialData = AddSubjectAndSession(trialData, subjectNum, sessionNum);

        % Remove duplicate rows from trial data.
        if height(trialData) > 1
            if isequal(trialData(height(trialData), :), trialData(height(trialData) - 1, :))
                trialData(height(trialData), :) = [];
            end
        end

        % Save the processed trial data to the current processed data path.
        SaveDataCheckCompleteFile('trialData', trialData, currentProcessedDataPath);
        % SaveDataCheckCompleteFile('abortedTrialData', abortedTrialData, processedDataPath);
        % SaveDataToMat(trialData, processedDataPath, 'TrialData');
    end
end

%% Process Block Data
function blockData = processBlockData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum, exptType, ignoreBlockCondition)
    % Attempt to load block data from a processed data file.
    blockData = LoadDataCheckCompleteFile('blockData', currentProcessedDataPath, forceProcessAllData);

    % If block data is not loaded, read it from raw data files.
    if isempty(blockData)
        % Read block data from raw data files.
        blockData = ReadDataFiles([runtimeDataPath filesep 'BlockData'], '*BlockData.txt', 'ImportOptions', {'delimiter', '\t', 'TreatAsEmpty',{'null'}});
        blockData = AddSubjectAndSession(blockData, subjectNum, sessionNum);

        % Save the processed block data to the current processed data path.
        SaveDataCheckCompleteFile('blockData', blockData, currentProcessedDataPath);
    end

    % Check for duplicate rows in block data and handle bad block information.
    if height(blockData) > 1
        if isequal(blockData(height(blockData), :), blockData(height(blockData) - 1, :))
            badBlock = blockData.Block(height(blockData));
        else
            badBlock = [];
        end
    else
        badBlock = [];
    end

    % Process block data further based on experiment type and configuration.
    if strcmp(exptType, 'FLU_GL') && (~ismember('MeanPositiveTokens', blockData.Properties.VariableNames))
        blockjson = jsondecode(fileread([runtimeDataPath filesep 'SessionSettings' filesep 'BlockDefs.json']));
        blockData = AddTokenInfo(blockData, blockjson);

        % Save the updated block data to the current processed data path.
        SaveDataCheckCompleteFile('blockData', blockData, currentProcessedDataPath);
    end

    % Process block data based on additional conditions if specified.
    if (~ismember('ID_ED', blockData.Properties.VariableNames) || forceProcessAllData) && ~ignoreBlockCondition
        blockData = AddBlockType(blockData, exptType, runtimeDataPath);
        % Save the updated block data to the current processed data path.
        SaveDataCheckCompleteFile('blockData', blockData, currentProcessedDataPath);
    end
end

%% Process Frame Data
function frameData = processFrameData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum)
    % Attempt to load frame data from a processed data file.
    frameData = LoadDataCheckCompleteFile('frameData', currentProcessedDataPath, forceProcessAllData);

   
    if exist([runtimeDataPath filesep 'FrameData'], 'dir')
        % Enter FrameData folder when parsing task FrameData
        path = [runtimeDataPath filesep 'FrameData'];
    else
        % Use only runtimeDataPath when parsing TaskSelectionData (there is only one FrameData file)
        path = runtimeDataPath;
    end
    
    % If frame data is not loaded or has only one column, read it from raw data files.
    if isempty(frameData) || width(frameData) == 1
        % Read frame data from raw data files.
        frameData = ReadDataFiles(path, '*FrameData*.txt', 'importOptions', {'delimiter', '\t', 'TreatAsEmpty',{'None'}});
        frameData = AddSubjectAndSession(frameData, subjectNum, sessionNum);

        % Save the processed frame data to the current processed data path.
        SaveDataCheckCompleteFile('frameData', frameData, currentProcessedDataPath);
    end
end

%% Process Gaze Data
function gazeData = processGazeData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum)
    % Display a message indicating that gaze data is being loaded.
    fprintf('\tLoading gaze data from file.\n');
    
    % Attempt to load gaze data from a processed data file.
    gazeData = LoadDataCheckCompleteFile('gazeData', currentProcessedDataPath, forceProcessAllData);

    if exist([runtimeDataPath filesep 'GazeData'], 'dir')
        % Enter GazeData folder when parsing Task GazeData
        path = [runtimeDataPath filesep 'GazeData'];
    else
        % Use only runtimeDataPath when parsing TaskSelectionData (there is only one GazeData file)
        path = runtimeDataPath;
    end

    % If gaze data is not loaded or has only one column, read it from raw data files.
    if isempty(gazeData) || width(gazeData) == 1
        % Read gaze data from raw data files.
        gazeData = ReadDataFiles(path, '*GazeData*.txt', 'importOptions', {'delimiter', '\t', 'TreatAsEmpty',{'None'}});
        gazeData = AddSubjectAndSession(gazeData, subjectNum, sessionNum);

        % Save the processed gaze data to the current processed data path.
        SaveDataCheckCompleteFile('gazeData', gazeData, currentProcessedDataPath);
    end
end

%% Process Serial Data
function [serialSentData, serialRecvData] = processSerialData(runtimeDataPath, forceProcessAllData, currentProcessedDataPath, subjectNum, sessionNum)
    % Load sent and received serial data if available.
    serialSentData.Raw = LoadDataCheckCompleteFile('serialSentData', currentProcessedDataPath, forceProcessAllData);
    serialRecvData.Raw = LoadDataCheckCompleteFile('serialRecvData', currentProcessedDataPath, forceProcessAllData);
   
    if exist([runtimeDataPath filesep 'SerialSentData'], 'dir')
        % Enter SerialSentData folder when parsing Task SerialSentData
        sentPath = [runtimeDataPath filesep 'SerialSentData'];
    else
        % Use only runtimeDataPath when parsing TaskSelectionData (there is only one SerialSentData file)
        sentPath = runtimeDataPath;
    end
    
    if exist([runtimeDataPath filesep 'SerialRecvData'], 'dir')
        % Enter SerialSentData folder when parsing Task SerialSentData
        recvPath = [runtimeDataPath filesep 'SerialRecvData'];
    else
        % Use only runtimeDataPath when parsing TaskSelectionData (there is only one SerialSentData file)
        recvPath = runtimeDataPath;
    end
    
    % If sent serial data is not loaded or has only one column, read it from raw data files.
    if isempty(serialSentData) || width(serialSentData) == 1
        % Read sent serial data from raw data files.
        serialSentData.Raw = ReadDataFiles(sentPath, '*SerialSentData*.txt', 'importOptions', {'delimiter', '\t', 'TreatAsEmpty',{'None'}});
        serialSentData.Raw = AddSubjectAndSession(serialSentData.Raw, subjectNum, sessionNum);
        
        % Save the processed sent serial data to the current processed data path.
        SaveDataCheckCompleteFile('serialSentData', serialSentData, currentProcessedDataPath);
    end
    
    % If received serial data is not loaded or has only one column, read it from raw data files.
    if isempty(serialRecvData) || width(serialRecvData) == 1
        % Read received serial data from raw data files.
        serialRecvData.Raw = ReadDataFiles(recvPath, '*SerialRecvData*.txt', 'importOptions', {'delimiter', '\t', 'TreatAsEmpty',{'None'}});
        serialRecvData.Raw = AddSubjectAndSession(serialRecvData.Raw, subjectNum, sessionNum);
        % serialRecvData = ParseSyncboxData(serialRecvData);
        serialRecvData.Analog = euUSE_parseSerialRecvDataAnalog(serialRecvData.Raw);
        [ serialRecvData.SynchA serialRecvData.SynchB serialRecvData.RwdA serialRecvData.RwdB serialRecvData.EventCodes] = euUSE_parseSerialRecvData( serialRecvData.Raw, 'word' );
        
        
        % Save the processed received serial data to the current processed data path.
        SaveDataCheckCompleteFile('serialRecvData', serialRecvData, currentProcessedDataPath);
    end
end

%% Helper Functions

function data = AddSubjectAndSession(data, subjectNum, sessionNum)
    if ~ismember('SessionNum', data.Properties.VariableNames)
        data = [table(repmat(sessionNum, height(data), 1), 'VariableNames', {'SessionNum'}) data];
    end
    if ~ismember('SubjectNum', data.Properties.VariableNames)
        data = [table(repmat(subjectNum, height(data), 1), 'VariableNames', {'SubjectNum'}) data];
    end
end

function SaveDataCheckCompleteFile(varname, data, dataFolderPath)
    SaveData(dataFolderPath, [upper(varname(1)) varname(2:end) '.mat'], varname, data);
end

function data = LoadDataCheckCompleteFile(varname, dataFolderPath, forceProcessAllData)
     global singleDataFileName singleDataFileVarNames multiDataFileVarNames; 

    iFile = strcmp(multiDataFileVarNames, [upper(varname(1)) varname(2:end) '.mat']);
    iVar = strcmp(singleDataFileVarNames, varname);

    if ~isempty(iFile) && sum(iFile) > 0 && ~forceProcessAllData
        filename = multiDataFileVarNames{iFile};
        data = LoadData(dataFolderPath, filename, varname);
    elseif ~isempty(iVar) && sum(iVar) > 0 && ~forceProcessAllData
        data = LoadData(dataFolderPath, singleDataFileName, varname);
    else
        data = [];
    end
end

function SaveData(dataFolderPath, dataFilename, varname, data)
    if ~exist(dataFolderPath, 'dir')
        mkdir(dataFolderPath);
    end
    fprintf(['\tSaving ' varname ' to file.\n']);
    S.(varname) = data;

    if ~exist([dataFolderPath filesep dataFilename], 'file')
        save([dataFolderPath filesep dataFilename], '-struct', 'S')
    else
        save([dataFolderPath filesep dataFilename], '-append', '-struct', 'S');
    end
end

function data = LoadData(dataFolderPath, dataFilename, varname)
    if ~exist([dataFolderPath filesep dataFilename], 'file')
        error(['Tried to load the variable ' varname ' from ' [dataFolderPath filesep dataFilename] ' but this file does not exist.']);
    else
        try
            fprintf(['\tLoading ' varname ' from file.\n']);
            S = load([dataFolderPath filesep dataFilename], varname);
            data = S.(varname);
        catch e
            error([e.identifier '/t' e.message]);
        end
    end
end

function SaveDataToMat(T, dataFolderPath, Tname)
    % Convert the table to an array
    data = table2struct(T);
    % Get the field names of the struct
    field_names = fieldnames(data);

    % Iterate through the fields of the struct
    for i = 1:numel(field_names)
        % Get the current field name
        field_name = field_names{i};

        % Get the current field value
        field_value = data.(field_name);

        % Save the field value to a .mat file with the same name as the field
        savepath = [dataFolderPath filesep Tname];
        if ~exist(savepath, 'dir')
            mkdir(savepath);
        end
        save([savepath filesep field_name], 'field_value');
    end
end

function CleanupFun(path)
    cd(path);
    warning('on', 'all');
end

function matches = findTaskNames(dataFolder)
    % Use the dir function to list all items in the directory
    dirInfo = dir(dataFolder);

    % Filter out only the directories (folders)
    folderNames = {dirInfo([dirInfo.isdir]).name};

    % Remove '.' (current directory) and '..' (parent directory) from the list
    folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

    % Regular expression pattern to match TaskXXXX_Name
    pattern = '^Task\d{4}_\w+';

    % Use cellfun and regexp to extract task names
    matches = regexp(folderNames, pattern, 'match');
    matches = vertcat(matches{:});  % Convert cell array of cell arrays to a single cell array

    % Display the extracted TaskXXXX_Name strings
    disp('Extracted Task Names:');
    disp(matches);
end
