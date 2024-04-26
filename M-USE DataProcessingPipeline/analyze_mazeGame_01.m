
% ---
% --- This script contains some basic analysis of the Maze task
% ---
%
% Copyright (c) 2023, Thilo Womelsdorf
%
% This file is part of the Multitask Unified Suite for Experiments (M-USE),
% see http://m-use.psy.vanderbilt.edu for documentation and details.
%
%    Matlab scripts associated with M-USE are free software:
%    you can redistribute it and/or modify it under the terms of the
%    GNU General Public License as published by the Free Software
%    Foundation,  either version 3 of the License, or (at your option)
%    any later version.
%
%    M-USE matlab scripts are distributed in the hope that it will be
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
%    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    A copy of the GNU General Public License can be found
%    at http://www.gnu.org/licenses/


do_LoadData = 1;
do_AnalysisBasic = 1;
do_AnalysisBasicPlot = 1;
do_AcrossSessionAnalysis_SortedByMazeType = 0;
do_AcrossSessionAnalysis_SortedBySession = 0;
do_AcrossSessionAnalysis_SortedByWeekday = 0;
do_AcrossSessionAnalysis_LearningCurve = 1;
do_AcrossSessionAnalysis_TurnVsStraightErrorRate = 0;
do_AcrossSessionAnalysis_TowardsVsAwayEndTile = 0;
do_AcrossSessionAnalysis_EarlyVsLate = 0;
do_AcrossSessionAnalysis_ReactionTimes = 0;

do_WithinSessionAnalysis = 0; % NEEDS TO BE UPDATED

DO_SAVEFIGURES = 0;

PLOT_INDIVIDUALSESSIONS = 0;

% MIN_TURNS = 4; % NEEDS TO BE UPDATED, SET TO 0 TO IGNORE TRIAL DATA FILTERING
% MIN_LENGTH = 14; % NEEDS TO BE UPDATED
% MAX_TURNS = 0;
% MAX_LENGTH = 0;

% --- --- --- --- --- --- --- ---
% --- folder with preprocecced data
% --- --- --- --- --- --- --- ---
HOME_FOLDER = ['C:\Users\Sorti\OneDrive\Documents\GitHub\M-USE_Analysis\M-USE DataProcessingPipeline' filesep];
%'/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS'
MUSEMATFOLDERNames = {}; iResultFolder = ''; iResultFile = '';

% SessionID_MZ = {'Frey_MZ_all_01'};
% SessionID_MZ = {'Wotan_MZ_all_01'};

   SessionID_MZ = {'Frey_MZ_all_02'};
   % SessionID_MZ = {'Wotan_MZ_all_02'};

% MUSEMATFOLDERNames = 'MUSEMAT01_VS_EC_FL_MZG_Frey'
switch SessionID_MZ{1}
    case {'Frey_MZ_all_01'}
        MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_FL_MZG_Frey'
        % MUSEMATFOLDERNames{2} = 'MUSEMAT01_VS_EC_FL_MRT_Frey'
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Frey_01_20240322'; %'WM_01_20231205';%'WM_01_20230907';
    case {'Frey_MZ_all_02'}
        MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_MG_Frey'
        % MUSEMATFOLDERNames{2} = 'MUSEMAT01_VS_EC_FL_MRT_Frey'
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Frey_01_20240402'; %'WM_01_20231205';%'WM_01_20230907';
    case {'Wotan_MZ_all_01'}
        MUSEMATFOLDERNames{1} =  'MUSEMAT01_WM_EC_FL_MZG_Wotan';
        %MUSEMATFOLDERNames{2} =  'MUSEMAT01_VS_EC_FL_MRT_Wotan';
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Wotan_01_20240322';%'WM_01_20230907';
    case {'Wotan_MZ_all_02'}
        MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_MG_Wotan';
        % MUSEMATFOLDERNames{2} = 'MUSEMAT01_VS_EC_FL_MRT_Frey'
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Wotan_01_20240402'; %'WM_01_20231205';%'WM_01_20230907';
end
iMonkey = SessionID_MZ{1}(1:findstr(SessionID_MZ{1},'_')-1);

% --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
% --- The files generates by this script and the folder name for them:
iResultFileMetrics = [iResultFile '_Metrics'];
RESULTFOLDER = [HOME_FOLDER filesep iResultFolder '_MAT']; if ~exist(RESULTFOLDER), mkdir(RESULTFOLDER),end
METRICSFOLDER  = [RESULTFOLDER '_METRICS']; if ~exist(METRICSFOLDER), mkdir(METRICSFOLDER),end
FIGURE_Folder = [RESULTFOLDER '_FIG'];    if ~exist(FIGURE_Folder), mkdir(FIGURE_Folder), end%
%
% % --- --- --- --- --- --- --- ---
% % --- folder with preprocecced data
% % --- --- --- --- --- --- --- ---
% MUSEMATFOLDERNames = {};
% % MUSEMATFOLDERNames = 'MUSEMAT01_VS_EC_FL_MZG_Frey'
% MUSEMATFOLDERNames{1} =  'MUSEMAT01_WM_EC_FL_MZG_Wotan';
% MUSEMATFOLDERNames{2} =  'MUSEMAT01_VS_EC_FL_MRT_Wotan';
%
% % --- The files generates by this script and the folder name for them:
% iResultFolder = 'RES01_WM_EC_FL_MZG_VS_Wotan';
% iResultFile   = 'MZ_01_20230907';
% iResultFileMetrics = [iResultFile '_Metrics'];
%
% % --- save collected results from do_LoadData
% RESULTFOLDER = [pwd filesep iResultFolder '_MAT']; if ~exist(RESULTFOLDER), mkdir(RESULTFOLDER),end
%
% % --- save analysis results and metrics
% METRICSFOLDER  = [RESULTFOLDER '_METRICS']; if ~exist(METRICSFOLDER), mkdir(METRICSFOLDER),end


% --- --- --- --- --- --- --- ---
% --- load and collect Data
% --- --- --- --- --- --- --- ---
if do_LoadData == 1

    % --- --- --- --- --- --- --- --- --- --- --- ---
    % --- collect all preprocessed data files (from those folders specified:
    % --- --- --- --- --- --- --- --- --- --- --- ---
    datasets = {};
    datasetsDIR = {};
    for iRF=1:length(MUSEMATFOLDERNames)
        DATAFOLDER = [pwd filesep MUSEMATFOLDERNames{iRF}];
        dirinfo = dir(DATAFOLDER);
        if isempty(dirinfo), sprintf('could not open %s\n',DATAFOLDER), return, end
        for j=1:length(dirinfo)
            if (dirinfo(j).isdir) | strcmp(dirinfo(j).name(1),'.' ),  continue, end
            if ~isempty(findstr(dirinfo(j).name,'DAT01_'))
                cL = length(datasets) + 1;
                datasetsDIR{cL} = DATAFOLDER;
                datasets{cL} = dirinfo(j).name(1:end-4);
            end
        end
    end

    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    % --- trialDataErrorAnalysis
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    %DO_TRIALDATA_LEVEL_ANALYSIS = 1; if DO_TRIALDATA_LEVEL_ANALYSIS

    % Initialize an empty struct array for results
    res = struct();

    for iD = 1:length(datasets)
        % Initialize dataset-specific information
        res(iD).dataset = datasets{iD};
        res(iD).datafolder = datasetsDIR{iD};
        res(iD).log = {};
        res(iD).data = struct(); % This will hold data for each label

        disp(sprintf('loading iD=%d, %s', iD, res(iD).dataset));
        in = load(fullfile(res(iD).datafolder, res(iD).dataset));
        disp('.');


        % Find indices of specific labels within the taskLabel field
        labels = {'Maze1_', 'MazeGuided_', 'Maze1Repeat_'};
        for label = labels
            idx_all = find(contains(in.dat.taskLabel, label));
            if isempty(idx_all)
                continue; % Skip to the next iteration of the loop
            end
            for idx = idx_all
                if idx > length(in.dat.trialData)
                    res(iD).log{end+1} = sprintf('Index exceeds number of trials... iD=%d, %s', iD, datasets{iD});
                    continue;
                end
                if isempty(idx)
                    res(iD).log{end+1} = sprintf('Index exceeds number of trials... iD=%d, %s', iD, datasets{iD});
                    continue;
                end
                % Check if there's data for the current index
                if ~isempty(in.dat.trialData{idx})
                    % Initialize the structure for the current label if it does not exist
                    labelStr = strrep(label{1}, '_', '');
                    if ~isfield(res(iD).data, labelStr)
                        res(iD).data.(labelStr) = struct();
                    end

                    if ~isempty(in.dat.frameData)
                        frameData = in.dat.frameData{idx};
                    end
                    trialData = in.dat.trialData{idx};
                    blockData = in.dat.blockData{idx};
                    blockDef = table2struct(in.dat.cfg_BlockDef{idx});
                    res(iD).subjectName = trialData(1).SubjectID;
                    iT = 0;
                    sessionDates = cellfun(@(x) regexp(x, 'Session_(\d+_\d+)', 'tokens'), {res(iD).dataset}, 'UniformOutput', false);
                    sessionDates = cellfun(@(x) x{1}{1}, sessionDates, 'UniformOutput', false); % Convert tokens to strings


                    % Determine the correct label based on whether the block ends with 'N'
                    newBlockIndices = find(arrayfun(@(x) endsWith(x.BlockName, 'N'), blockDef));
                    repeatBlockIndices = find(arrayfun(@(x) endsWith(x.BlockName, 'R'), blockDef));
                    currentSessionData = struct(); % data for the current session

                    currentSessionData.sessionDatesFormatted = datetime(sessionDates, 'InputFormat', 'MM_dd', 'Format', 'MM/dd'); % Convert to datetime
                    res(iD).sessionDatesFormatted = currentSessionData.sessionDatesFormatted;

                    if(~isempty(newBlockIndices))
                        newTrialData = trialData(arrayfun(@(x) ismember(x.BlockCount, newBlockIndices), trialData));
                        currentSessionData = AppendCurrentSessionData(newTrialData, currentSessionData, blockDef);
                        res(iD).data.Maze2 = currentSessionData;
                    end
                    if(~isempty(repeatBlockIndices))
                        repeatTrialData = trialData(arrayfun(@(x) ismember(x.BlockCount, repeatBlockIndices), trialData));
                        currentSessionData = AppendCurrentSessionData(repeatTrialData, currentSessionData, blockDef);
                        res(iD).data.Maze1Repeat = currentSessionData;
                    end
                    if(isempty(newBlockIndices) && isempty(repeatBlockIndices))
                        currentSessionData = AppendCurrentSessionData(trialData, currentSessionData, blockDef);
                        res(iD).data.Maze1 = currentSessionData;
                    end


                else
                    res(iD).log{end+1} = sprintf('Trial data is empty for this label... iD=%d, %s', iD, datasets{iD});
                end
            end
        end
    end

    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    % --- save analysis results in RES file
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    save(fullfile(RESULTFOLDER, iResultFile), 'res', '-v7.3');
    fprintf('saved analysis of n=%d datasets in %s.', iD, iResultFile);
end


% - #Time (in sec.) to reach target-Tile (From first touch to target-Tile touch)
% - #Total Number of touches to reach goal
% - #Normalized number of touches to end goal () Ttal touches needed divded by optimal number required)
% - #Total Number of rule-breaking errors
% - #Proportion of rule-breaking errors
% - #Total Number of rule-abiding errors
% - #Proportion of rule-abiding errors



% --- --- --- --- --- --- --- --- --- --- --- --- ---
% --- do some basic analysis
%% --- --- --- --- --- --- --- --- --- --- --- --- ---
if do_AnalysisBasic == 1

    % Load results if they are not yet loaded
    if ~exist('res', 'var')
        load(fullfile(RESULTFOLDER, iResultFile))
    end
    if do_AcrossSessionAnalysis_SortedByMazeType == 1 % compile all unique combinations of turns and length across sessions
        allMazes = [];
        for i = 1:length(res)
            if isfield(res(i).data, 'Maze1') && ~isempty(fieldnames(res(i).data.Maze1))
                allMazes = [allMazes; res(i).data.Maze1.mazeTurnsLength];
            end
            if isfield(res(i).data, 'Maze1Repeat') && ~isempty(fieldnames(res(i).data.Maze1Repeat))
                allMazes = [allMazes; res(i).data.Maze1Repeat.mazeTurnsLength];
            end
            if isfield(res(i).data, 'Maze2') && ~isempty(fieldnames(res(i).data.Maze2))
                allMazes = [allMazes; res(i).data.Maze2.mazeTurnsLength];
            end
        end

        uniqueMazeConfigurations = unique(allMazes, 'rows');
    end

    if (do_AcrossSessionAnalysis_LearningCurve == 1 || do_WithinSessionAnalysis == 1 || do_AcrossSessionAnalysis_EarlyVsLate || do_AcrossSessionAnalysis_TurnVsStraightErrorRate || do_AcrossSessionAnalysis_TowardsVsAwayEndTile || do_AcrossSessionAnalysis_ReactionTimes)
        allTrials= [];
        for i = 1:length(res)
            if isfield(res(i).data, 'Maze1') && ~isempty(fieldnames(res(i).data.Maze1))
                allTrials = [allTrials; res(i).data.Maze1.trialInBlock'];
            end
            if isfield(res(i).data, 'Maze1Repeat') && ~isempty(fieldnames(res(i).data.Maze1Repeat))
                allTrials = [allTrials; res(i).data.Maze1Repeat.trialInBlock'];
            end
            if isfield(res(i).data, 'Maze2') && ~isempty(fieldnames(res(i).data.Maze2))
                allTrials = [allTrials; res(i).data.Maze2.trialInBlock'];
            end
        end

        uniqueTrialCounts = unique(allTrials, 'rows');
    end


    % Initialize the empty arrays to store data for the respective analysis
    % type
    mazeTypes = {'Maze1', 'Maze2', 'Maze1Repeat', 'Aggregate'};
    standardMetrics = {'normalizedTotalErrors', 'normalizedMazeDurations', ...
        'normalizedRuleBreakingErrors', 'normalizedRuleAbidingErrors', ...
        'normalizedPerseverativeErrors'};

    if do_AcrossSessionAnalysis_SortedByMazeType == 1
        metrics = [standardMetrics, {'numTurns', 'numLength'}];
        metricsLength = size(uniqueMazeConfigurations, 1);
    elseif do_AcrossSessionAnalysis_SortedBySession == 1
        metrics = [standardMetrics, {'sessionDate'}];
        metricsLength = length(res);
    elseif do_AcrossSessionAnalysis_SortedByWeekday == 1
        metrics = [standardMetrics, {'dayOfTheWeek'}];
        metricsLength = 5;
    elseif (do_AcrossSessionAnalysis_LearningCurve == 1 || do_WithinSessionAnalysis == 1 || do_AcrossSessionAnalysis_EarlyVsLate)
        metrics = [standardMetrics, {'trialInBlock'}];
        metricsLength = size(uniqueTrialCounts, 1);
    elseif(do_AcrossSessionAnalysis_TurnVsStraightErrorRate == 1 || do_AcrossSessionAnalysis_TowardsVsAwayEndTile)
        metrics = {'ruleBreakingErrors', 'ruleAbidingErrors', ...
            'perseverativeErrors', 'retouchErrors', 'backTrackErrors'};
    elseif (do_AcrossSessionAnalysis_ReactionTimes == 1)
        metrics = {'correctToCorrect', 'correctToError', ...
            'errorToRetouch', 'retouchToCorrect', 'retouchToError', 'errorToError'};
        metricsLength = 1;
    end


    % Initialize the structure
    metrics_mz = struct();

    % Initialize the structure for 'turn' and 'straight' with NaN arrays
    for i = 1:length(mazeTypes)
        mazeType = mazeTypes{i};
        for j = 1:length(metrics)
            metric = metrics{j};
            if (do_AcrossSessionAnalysis_TurnVsStraightErrorRate == 1 || do_AcrossSessionAnalysis_TowardsVsAwayEndTile == 1)
                if do_AcrossSessionAnalysis_TurnVsStraightErrorRate == 1
                    % Pre-allocate with NaNs, one row for each unique trial count
                    metrics_mz.(mazeType).turn.(metric) = cell(size(uniqueTrialCounts, 1), 1);
                    metrics_mz.(mazeType).straight.(metric) = cell(size(uniqueTrialCounts, 1), 1);
                end
                if do_AcrossSessionAnalysis_TowardsVsAwayEndTile == 1
                    metrics_mz.(mazeType).towards.(metric) = cell(size(uniqueTrialCounts, 1), 1);
                    metrics_mz.(mazeType).away.(metric) = cell(size(uniqueTrialCounts, 1), 1);
                end
            else
                % Other initializations for different types of analysis
                metrics_mz.(mazeType).(metric) = cell(metricsLength, 1);
            end
        end
    end


    if(do_WithinSessionAnalysis ==1 )
        SESSION_DATE = "03/18";
        sessionIdx = find(res.sessionDatesFormatted == SESSION_DATE);

        maze1Data = res(sessionIdx).data.Maze1;
        repeatMaze1Data = res(sessionIdx).data.Maze1Repeat;
        metrics_mz = AppendMetrics(metrics_mz, 'Maze1', maze1Data, iL, filteredMaze1Indices);
        metrics_mz.Maze1.sessionDate{sessionIdx} = maze1Data.sessionDatesFormatted;

        metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMaze1Data, iL, filteredMaze1Indices);

        metrics_mz.Maze1Repeat.sessionDate{sessionIdx} = repeatMaze1Data.sessionDatesFormatted;

    else
        for sessionIdx = 1:length(res)
            maze1Data = res(sessionIdx).data.Maze1;
            if isfield(res(sessionIdx).data, 'Maze1Repeat')
                repeatMaze1Data = res(sessionIdx).data.Maze1Repeat;
            end
            if isfield(res(sessionIdx).data, 'Maze2')
                maze2Data = res(sessionIdx).data.Maze2;
            end

            if do_AcrossSessionAnalysis_SortedByMazeType == 1
                for iL = 1:size(uniqueMazeConfigurations, 1)
                    filteredMaze1Indices = find(maze1Data.mazeTurnsLength(:,1) ==  uniqueMazeConfigurations(iL,1) & maze1Data.mazeTurnsLength(:,2) == uniqueMazeConfigurations(iL,2));
                    % if ~isempty(filteredMaze1Indices)
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1', maze1Data, iL, filteredMaze1Indices);
                    metrics_mz.Maze1.numTurns{iL} = uniqueMazeConfigurations(iL,1);
                    metrics_mz.Maze1.numLength{iL} = uniqueMazeConfigurations(iL,2);

                    filteredMaze1RepeatIndices = find(repeatMaze1Data.mazeTurnsLength(:,1) ==  uniqueMazeConfigurations(iL,1) & repeatMaze1Data.mazeTurnsLength(:,2) == uniqueMazeConfigurations(iL,2));
                    % if ~isempty(filteredMaze1RepeatIndices)
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMaze1Data, iL, filteredMaze1RepeatIndices);
                    metrics_mz.Maze1Repeat.numTurns{iL} = uniqueMazeConfigurations(iL,1);
                    metrics_mz.Maze1Repeat.numLength{iL} = uniqueMazeConfigurations(iL,2);

                    filteredMaze2Indices = find(maze2Data.mazeTurnsLength(:,1) ==  uniqueMazeConfigurations(iL,1) & maze2Data.mazeTurnsLength(:,2) == uniqueMazeConfigurations(iL,2));
                    % if ~isempty(filteredMaze1RepeatIndices)
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze2', maze2Data, iL, filteredMaze2Indices);
                    metrics_mz.Maze2.numTurns{iL} = uniqueMazeConfigurations(iL,1);
                    metrics_mz.Maze2.numLength{iL} = uniqueMazeConfigurations(iL,2);
                end

            elseif do_AcrossSessionAnalysis_SortedBySession == 1
                metrics_mz = AppendMetrics(metrics_mz, 'Maze1', maze1Data, sessionIdx, 1:height(maze1Data.mazeTurnsLength));
                metrics_mz.Maze1.sessionDate{sessionIdx} = maze1Data.sessionDatesFormatted;

                metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMaze1Data, sessionIdx, 1:height(repeatMaze1Data.mazeTurnsLength));
                metrics_mz.Maze1Repeat.sessionDate{sessionIdx} = repeatMaze1Data.sessionDatesFormatted;

                metrics_mz = AppendMetrics(metrics_mz, 'Maze2', maze2Data, sessionIdx, 1:height(maze2Data.mazeTurnsLength));
                metrics_mz.Maze2.sessionDate{sessionIdx} = maze2Data.sessionDatesFormatted;

            elseif do_AcrossSessionAnalysis_SortedByWeekday == 1
                daysOfTheWeek = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
                for iDay = 1:5
                    day = daysOfTheWeek(iDay);
                    filteredMaze1Indices = find(strcmpi(maze1Data.dayOfTheWeek(1,:), day));
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1', maze1Data, iDay, filteredMaze1Indices);
                    metrics_mz.Maze1.dayOfTheWeek(iDay) = day;

                    filteredMaze1RepeatIndices = find(strcmpi(repeatMaze1Data.dayOfTheWeek(1,:), day));
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMaze1Data, iDay, filteredMaze1RepeatIndices);

                    filteredMaze2Indices = find(strcmpi(maze2Data.dayOfTheWeek(1,:), day));
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze2', maze2Data, iDay, filteredMaze2Indices);

                end
            elseif (do_AcrossSessionAnalysis_LearningCurve == 1 || do_AcrossSessionAnalysis_EarlyVsLate == 1)

                for iTrialInBlock = 1:height(uniqueTrialCounts)
                    filteredMaze1Indices = find(maze1Data.trialInBlock == iTrialInBlock);
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1', maze1Data, iTrialInBlock, filteredMaze1Indices);

                    if isfield(res(sessionIdx).data, 'Maze1Repeat') && ~isempty(fieldnames(res(sessionIdx).data.Maze1Repeat))
                        filteredMaze1RepeatIndices = find(repeatMaze1Data.trialInBlock == iTrialInBlock);
                        metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMaze1Data, iTrialInBlock, filteredMaze1RepeatIndices);
                    end

                    if isfield(res(sessionIdx).data, 'Maze2')
                        filteredMaze2Indices = find(maze2Data.trialInBlock == iTrialInBlock);
                        metrics_mz = AppendMetrics(metrics_mz, 'Maze2', maze2Data, iTrialInBlock, filteredMaze2Indices);
                    end
                end
                metrics_mz.Maze1.trialInBlock = num2cell(1:uniqueTrialCounts(end));
                metrics_mz.Maze1Repeat.trialInBlock = num2cell(1:uniqueTrialCounts(end));
                metrics_mz.Maze2.trialInBlock = num2cell(1:uniqueTrialCounts(end));


            elseif (do_AcrossSessionAnalysis_TurnVsStraightErrorRate == 1 || do_AcrossSessionAnalysis_TowardsVsAwayEndTile == 1)
                % Process each maze data set
                metrics_mz = ProcessSpatialMazeData(metrics_mz, maze1Data, 'Maze1', do_AcrossSessionAnalysis_TurnVsStraightErrorRate, do_AcrossSessionAnalysis_TowardsVsAwayEndTile);
                metrics_mz = ProcessSpatialMazeData(metrics_mz, repeatMaze1Data, 'Maze1Repeat', do_AcrossSessionAnalysis_TurnVsStraightErrorRate, do_AcrossSessionAnalysis_TowardsVsAwayEndTile);
                metrics_mz = ProcessSpatialMazeData(metrics_mz, maze2Data, 'Maze2', do_AcrossSessionAnalysis_TurnVsStraightErrorRate,do_AcrossSessionAnalysis_TowardsVsAwayEndTile);


            elseif ( do_AcrossSessionAnalysis_ReactionTimes == 1)
                if(~isfield(maze1Data, 'reactionTimes'))
                    continue;
                else
                    metrics_mz = ProcessReactionTimeData(metrics_mz, maze1Data, 'Maze1');
                    metrics_mz = ProcessReactionTimeData(metrics_mz, repeatMaze1Data, 'Maze1Repeat');
                    metrics_mz = ProcessReactionTimeData(metrics_mz, maze2Data, 'Maze2');
                end

            end


        end

    end

    metrics_mz = AggregateMetrics(metrics_mz);

end

% % % 
if do_AnalysisBasicPlot == 1

    % --- --- --- --- --- --- --- --- --- --- ---
    % --- load metrics if they are not yet loaded
    % --- --- --- --- --- --- --- --- --- --- ---

    if ~exist('res')
        load([RESULTFOLDER filesep iResultFileMetrics])
    end
    % disp('wait here'), return

    %  --- --- --- --- --- --- ---
    %  --- which mazes to look at ?
    %  --- --- --- --- --- --- ---


    if do_AcrossSessionAnalysis_SortedByMazeType == 1
        plot_AcrossSessionAnalysis_SortedByMazeType(metrics_mz, res);
    end

    if do_AcrossSessionAnalysis_SortedBySession == 1
        plot_AcrossSessionAnalysis_SortedBySession(metrics_mz, res);
    end

    if do_AcrossSessionAnalysis_SortedByWeekday == 1
        plot_AcrossSessionAnalysis_SortedByWeekday(metrics_mz, res);
    end

    if do_AcrossSessionAnalysis_LearningCurve == 1
        plot_AcrossSessionAnalysis_LearningCurve(metrics_mz, res, uniqueTrialCounts);
    end

    if(do_AcrossSessionAnalysis_TurnVsStraightErrorRate == 1)
        plot_AcrossSessionAnalysis_TurnVsStraightErrorRate(metrics_mz, res);
    end

    if(do_AcrossSessionAnalysis_TowardsVsAwayEndTile == 1)
        plot_AcrossSessionAnalysis_TowardsVsAway(metrics_mz, res);
    end

    if(do_AcrossSessionAnalysis_EarlyVsLate == 1)
        plot_AcrossSessionAnalysis_EarlyVsLate(metrics_mz, res);
    end

    if(do_AcrossSessionAnalysis_ReactionTimes == 1)
        plot_AcrossSessionAnalysis_ReactionTimes(metrics_mz, res);
    end


end



if(PLOT_INDIVIDUALSESSIONS)
    % Prepare a single figure with two subplots
    figure;
    % Extract relevant fields from both structs
    blockNumsInitial = metrics_mz1.blockNum;
    trialInBlockInitial = metrics_mz1.trialInBlock;
    totalErrorsInitial = metrics_mz1.totalErrors;
    normalizeDurationsInitial = metrics_mz1.normalizedDurati
    mazeDefNameInitial = metrics_mz1.mazeDefName;

    blockNumsRepeat = metrics_mz1Repeat.blockNum;
    trialInBlockRepeat = metrics_mz1Repeat.trialInBlock;
    totalErrorsRepeat = metrics_mz1Repeat.totalErrors;
    mazeDefNameRepeat = metrics_mz1Repeat.mazeDefName;
    % Define colors for each block
    colors = lines(max(length(unique(blockNumsInitial)), length(unique(blockNumsRepeat))));
    legendEntries = {}; % We will collect legend info in two steps now

    % Initialize arrays to hold line objects for legend
    hLines = [];

    for i = 1:length(unique(blockNumsInitial))
        blockIndicesInitial = find(blockNumsInitial == i);
        matchingBlockIndex = find(strcmp(mazeDefNameInitial{blockIndicesInitial(1)}, mazeDefNameRepeat), 1);

        if ~isempty(matchingBlockIndex)
            blockIndicesRepeat = find(blockNumsRepeat == blockNumsRepeat(matchingBlockIndex));

            nTiles = metrics_mz1.mazeTurnsLength(blockIndicesInitial(1), 2);  % Assuming consistent maze length
            [sortedTrialsInitial, sortIndexInitial] = sort(trialInBlockInitial(blockIndicesInitial));
            [sortedTrialsRepeat, sortIndexRepeat] = sort(trialInBlockRepeat(blockIndicesRepeat));

            % Calculations for normalized errors and durations (assuming normalization is intended)
            normalizedErrorsInitial = totalErrorsInitial(blockIndicesInitial(sortIndexInitial)) ./ nTiles;
            normalizedErrorsRepeat = totalErrorsRepeat(blockIndicesRepeat(sortIndexRepeat)) ./ nTiles;
            normalizedDurationsInitial = metrics_mz1.mazeDuration(blockIndicesInitial(sortIndexInitial)) ./ nTiles;
            normalizedDurationsRepeat = metrics_mz1Repeat.mazeDuration(blockIndicesRepeat(sortIndexRepeat)) ./ nTiles;

            % First subplot for normalized total errors
            subplot(1, 2, 1);
            hold on;
            scatter(sortedTrialsInitial, normalizedErrorsInitial, 'filled', 'CData', colors(i,:));
            hLine1 = plot(sortedTrialsInitial, normalizedErrorsInitial, 'Color', colors(i,:), 'LineWidth', 1.5);
            scatter(sortedTrialsRepeat, normalizedErrorsRepeat, 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', 'none');
            hLine2 = plot(sortedTrialsRepeat, normalizedErrorsRepeat, '--', 'Color', colors(i,:), 'LineWidth', 1.5);

            % Add line handles and their descriptions to arrays for the legend
            hLines = [hLines, hLine1, hLine2];
            legendEntries{end+1} = ['Block ' num2str(i) ' Initial'];
            legendEntries{end+1} = ['Block ' num2str(i) ' Repeat'];


            % Second subplot for normalized maze durations
            subplot(1, 2, 2);
            hold on;
            scatter(sortedTrialsInitial, normalizedDurationsInitial, 'filled', 'CData', colors(i,:));
            plot(sortedTrialsInitial, normalizedDurationsInitial, 'Color', colors(i,:), 'LineWidth', 1.5);
            scatter(sortedTrialsRepeat, normalizedDurationsRepeat, 'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', 'none');
            plot(sortedTrialsRepeat, normalizedDurationsRepeat, '--', 'Color', colors(i,:), 'LineWidth', 1.5);
        end
    end

    % Finalize the first subplot with legend for lines only
    subplot(1, 2, 1);
    xlabel('Trial Number Within Block');
    ylabel('Normalized Total Errors');
    title('Normalized Total Errors');
    legend(hLines, legendEntries, 'Location', 'best');
    grid on;

    % Finalize the second subplot, repeating legend setup if necessary
    subplot(1, 2, 2);
    xlabel('Trial Number Within Block');
    ylabel('Normalized Total Errors');
    title('Normalized Total Errors');
    legend(hLines, legendEntries, 'Location', 'best');
    grid on;

    set(gcf, 'Position', [100, 100, 1400, 600]);
end



% Optionally save the figures
if DO_SAVEFIGURES
    currentDateNumber = datestr(now, 'yyyymmdd');
    figurefilename = sprintf('fig_mz_AcrossSessions_type1_%s.pdf', currentDateNumber);
    saveas(gcf, fullfile(FIGURE_Folder, figurefilename), 'pdf');
    close(gcf)
end

% %  --- Percent errors per block
%     iLegend = {'M1','M2','M3','M4','M5','M6'}
%     Z = resBlock.trlInBlockPerformance(:,1:5,2) ./ resBlock.trlInBlockPerformance(:,1:5,1);
%     X = 1:5;  % repmat(X,size(Z,1),1)
%
%     figure('Color','w'), hold on,
%     for iZ=1:size(Z,1)
%         plot(X+randi([-100 100],1,1)*0.001, Z(iZ,:),'o-','Color',iCol{iZ}), hold on
%     end
%     set(gca,'xtick',X,'tickdir','out')
%     ylabel('prop errors in block'), xlabel('trial in block')
%     legend(iLegend)



disp('done, return'), return


function metrics_mz = AppendMetrics(metrics_mz, mazeType, mazeData, metricsIndex, mazeDataIndex)

metrics_mz.(mazeType).normalizedTotalErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedTotalErrors{metricsIndex}, mazeData.totalErrors(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];
metrics_mz.(mazeType).normalizedMazeDurations{metricsIndex} = [metrics_mz.(mazeType).normalizedMazeDurations{metricsIndex}, mazeData.mazeDuration(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];
metrics_mz.(mazeType).normalizedRuleBreakingErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedRuleBreakingErrors{metricsIndex}, mazeData.ruleBreakingErrors(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];
metrics_mz.(mazeType).normalizedRuleAbidingErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedRuleAbidingErrors{metricsIndex}, mazeData.ruleAbidingErrors(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];

sumOfPerseverativeErrors = mazeData.perseverativeRuleBreakingErrors(mazeDataIndex) + ...
    mazeData.perseverativeRuleAbidingErrors(mazeDataIndex) + ...
    mazeData.perseverativeRetouchErroneous(mazeDataIndex);

metrics_mz.(mazeType).normalizedPerseverativeErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedPerseverativeErrors{metricsIndex}, sumOfPerseverativeErrors ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];

end

function metrics_mz = AggregateMetrics(metrics_mz)
% Iterate through each field name in Maze1
fields = fieldnames(metrics_mz.Maze1);

for f = 1:length(fields)
    fieldName = fields{f};
    % Initialize a cell to hold the aggregated data
    aggregatedData = {};

    % Aggregate data from Maze1, Maze1Repeat, and optionally Maze2
    % Check if the field is a struct (has subfields)
    if isstruct(metrics_mz.Maze1.(fieldName))
        % Iterate through subfields if present
        subfields = fieldnames(metrics_mz.Maze1.(fieldName));
        for sf = 1:length(subfields)
            subfieldName = subfields{sf};
            % Aggregate subfield data from Maze1, Maze1Repeat, and optionally Maze2
            aggregatedData = aggregateField(metrics_mz, fieldName, subfieldName);
            % Assign the aggregated data for the subfield
            metrics_mz.Aggregate.(fieldName).(subfieldName) = aggregatedData;
        end
    else
        % Directly aggregate data for fields without subfields
        aggregatedData = aggregateField(metrics_mz, fieldName);
        % Assign the aggregated data
        metrics_mz.Aggregate.(fieldName) = aggregatedData;
    end
end
end

% Helper function to aggregate data
function aggregatedData = aggregateField(metrics_mz, fieldName, subfieldName)
aggregatedData = {};
mazes = {'Maze1', 'Maze1Repeat', 'Maze2'};
data = {};

% Collect data from each Maze, if available
for m = 1:length(mazes)
    mazeName = mazes{m};
    if isfield(metrics_mz, mazeName)
        if nargin == 3 % Subfield present
            if isfield(metrics_mz.(mazeName).(fieldName), subfieldName)
                data{end+1} = metrics_mz.(mazeName).(fieldName).(subfieldName);
            end
        else
            if isfield(metrics_mz.(mazeName), fieldName)
                data{end+1} = metrics_mz.(mazeName).(fieldName);
            end
        end
    end
end

if ~isempty(data) && all(cellfun(@(x) isequal(size(x), size(data{1})), data))
    if iscell(data{1})
        % Number of rows in any of the cell arrays
        numRows = size(data{1}, 1);

        % Pre-allocate the aggregatedData cell array
        aggregatedData = cell(numRows, 1);

        % Loop over each row
        for i = 1:numRows
            % Extract the i-th element from each cell array and concatenate them
            aggregatedRow = cellfun(@(x) x{i}, data, 'UniformOutput', false);
            aggregatedData{i} = horzcat(aggregatedRow{:});
        end
    else
        % For numerical arrays, concatenate directly
        aggregatedData = horzcat(data{:});
    end
else
    numCols = size(data, 2);

    % Loop over each row
    for i = 1:numCols
        aggregatedData = [aggregatedData, data{i}];
    end

end

if(isempty(aggregatedData))
    error('Data size mismatch or empty data.');

end
end


%
% function metrics_mz = AggregateMetrics(metrics_mz)
% fields = fieldnames(metrics_mz.Maze1);
% subfield = {};
%
% for f = 1:length(fields)
%     fieldName = fields{f};
%     if isstruct(metrics_mz.Maze1.(fieldName))
%         subfield = fieldnames(metrics_mz.Maze1.(fieldName));
%         for sf = 1:length(subfield)
%             subfieldName = subfield{sf};
%             metrics_mz.Aggregate.(fieldName).(subfieldName) = [metrics_mz.Maze1.(fieldName).(subfieldName), metrics_mz.Maze1Repeat.(fieldName).(subfieldName)];
%         end
%     else
%
%
% metrics_mz.Aggregate.(fieldName) = [metrics_mz.Maze1.(fieldName), metrics_mz.Maze1Repeat.(fieldName)];
%     end
% end
% end

function adjacent = isAdjacent(tile1, tile2)
% Helper function to convert chess coordinate to Cartesian coordinate

[x1, y1] = getCartesian(tile1);
[x2, y2] = getCartesian(tile2);

dx = abs(x1 - x2);
dy = abs(y1 - y2);

adjacent = (dx == 1 && dy == 0) || (dx == 0 && dy == 1);
end

function [errorTypes, tileSpecificErrorCounts] = determineSelectionError(selectedTiles, mPath)
% Initialize
errorTypes = strings(1, length(selectedTiles)); % Use string array for error types
lastCorrectIndex = 1; % Last correct tile index in the path
consecutiveErrors = 0; % Count of consecutive errors
startedMaze = false; % Tracks if the start tile was correctly touched first
lastErrorTile = ""; % Last tile that caused an error

tileSpecificErrorCounts = struct("ruleBreakingErrors", zeros(6,6),"ruleAbidingErrors", zeros(6,6),"backTrackErrors", zeros(6,6),"retouchErrors", zeros(6,6), "perseverativeErrors", zeros(6,6));

% Initial target is the start tile
[targetX, targetY] = getCartesian(mPath{1});


% Iterate over selected tiles
for i = 1:length(selectedTiles)
    tile = selectedTiles{i};
    backTrackError = false;
    retouchCurrentTilePositionError = false;
    retouchCurrentTilePositionCorrect = false;
    ruleAbidingError = false;
    ruleBreakingError = false;
    correctNextTileChoice = false;

    % Check start condition
    if ~startedMaze
        if strcmp(tile, mPath{1})
            startedMaze = true;
            lastCorrectIndex = 1;
            correctNextTileChoice = true;
        else
            ruleBreakingError = true;
            consecutiveErrors = consecutiveErrors + 1;
        end
    else
        if strcmp(tile, mPath{min(lastCorrectIndex + 1, length(mPath))}) && consecutiveErrors == 0
            % Correct selection
            lastCorrectIndex = min(lastCorrectIndex + 1, length(mPath)); % sets the newly selected tile as the last correct idx
            correctNextTileChoice = true;
            consecutiveErrors = 0; % Reset errors
            lastErrorTile = ""; % Reset last error tile
        elseif strcmp(tile, mPath{lastCorrectIndex})
            % Retouching the last correct tile
            if consecutiveErrors > 0
                retouchCurrentTilePositionCorrect = true;
                consecutiveErrors = 0; % Reset errors

            else
                retouchCurrentTilePositionError = true;
            end
        elseif isAdjacent(tile, mPath{lastCorrectIndex}) && ~any(strcmp(mPath(1:lastCorrectIndex), tile))
            if(consecutiveErrors > 0 )
                % Rule-breaking error: failed to return to last correct
                % after an error
                ruleBreakingError = true;
                consecutiveErrors = consecutiveErrors + 1;
            else
                % Rule-abiding error: adjacent but incorrect
                ruleAbidingError = true;
                consecutiveErrors = consecutiveErrors + 1;
            end
        elseif any(strcmp(mPath(1:lastCorrectIndex), tile))
            backTrackError = true;
            consecutiveErrors = consecutiveErrors + 1;
        else
            % Rule-breaking error: not adjacent or repeating a previous error
            ruleBreakingError = true;
            consecutiveErrors = consecutiveErrors + 1;
        end
    end

    % Determine error type
    if ~isempty(lastErrorTile) && strcmp(lastErrorTile, tile) && ~correctNextTileChoice
        % Perseverative errors
        if backTrackError
            errorTypes(i) = "perseverativeBackTrackError";
            tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX) = tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX) + 1;
            tileSpecificErrorCounts.backTrackErrors(7 - targetY, targetX) = tileSpecificErrorCounts.backTrackErrors(7 - targetY, targetX) + 1;
            tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) = tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) + 1;
        elseif retouchCurrentTilePositionError
            errorTypes(i) = "perseverativeRetouchError";
            tileSpecificErrorCounts.retouchErrors(7 - targetY, targetX)  = tileSpecificErrorCounts.retouchErrors(7 - targetY, targetX) + 1;
            tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) = tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) + 1;

        elseif ruleAbidingError
            errorTypes(i) = "perseverativeRuleAbidingError";
            tileSpecificErrorCounts.ruleAbidingErrors(7 - targetY, targetX) = tileSpecificErrorCounts.ruleAbidingErrors(7 - targetY, targetX)  + 1;
            tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) = tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) + 1;

        else
            errorTypes(i) = "perseverativeRuleBreakingError";
            tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX)  = tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX)  + 1;
            tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) = tileSpecificErrorCounts.perseverativeErrors(7 - targetY, targetX) + 1;
        end
    else
        if retouchCurrentTilePositionCorrect
            errorTypes(i) = "retouchCorrect";
            [targetX, targetY] = getCartesian(mPath{lastCorrectIndex + 1});

        elseif correctNextTileChoice
            errorTypes(i) = "correct";
            [targetX, targetY] = getCartesian(mPath{min(lastCorrectIndex + 1, length(mPath))});

        elseif backTrackError
            errorTypes(i) = "backTrackError";
            tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX) = tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX)  + 1;
            tileSpecificErrorCounts.backTrackErrors(7 - targetY, targetX) = tileSpecificErrorCounts.backTrackErrors(7 - targetY, targetX) + 1;
            [targetX, targetY] = getCartesian(mPath{lastCorrectIndex});

        elseif retouchCurrentTilePositionError
            errorTypes(i) = "retouchError";
            tileSpecificErrorCounts.retouchErrors(7 - targetY, targetX) = tileSpecificErrorCounts.retouchErrors(7 - targetY, targetX)  + 1;
        elseif ruleAbidingError
            errorTypes(i) = "ruleAbidingError";
            tileSpecificErrorCounts.ruleAbidingErrors(7 - targetY, targetX) = tileSpecificErrorCounts.ruleAbidingErrors(7 - targetY, targetX) + 1;
            [targetX, targetY] = getCartesian(mPath{lastCorrectIndex});

        elseif ruleBreakingError
            errorTypes(i) = "ruleBreakingError";
            tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX) = tileSpecificErrorCounts.ruleBreakingErrors(7 - targetY, targetX)  + 1;
            [targetX, targetY] = getCartesian(mPath{lastCorrectIndex});

        end
    end

    % Update lastErrorTile if there was an error
    if ruleBreakingError || ruleAbidingError || backTrackError || retouchCurrentTilePositionError
        lastErrorTile = tile;
    end

    % Adjust for start condition
    if ~startedMaze
        errorTypes(i) = "ruleBreakingError"; % Overwrite for failure to start correctly
    end
end
end

function turns = getTurnsAlongPath(mPath)
turns = false(1, length(mPath));
for i = 2:length(mPath)-1
    turns(i) = checkIfTurn(mPath{i-1}, mPath{i}, mPath{i+1});
end
end
function isTowards = isMovingTowardsEnd(currentTile, nextTile, endTile)
% Convert chess coordinates to Cartesian
[currentX, currentY] = getCartesian(currentTile);
[nextX, nextY] = getCartesian(nextTile);
[endX, endY] = getCartesian(endTile);

% Vector from current to next tile
vectorCurrentToNext = [nextX - currentX, nextY - currentY];

% Vector from current to end tile
vectorCurrentToEnd = [endX - currentX, endY - currentY];

% Dot product between the two vectors
dotProduct = dot(vectorCurrentToNext, vectorCurrentToEnd);

% If dot product is positive, the direction is towards the end tile
isTowards = dotProduct > 0;
end

function isTurn = checkIfTurn(tile1, tile2, tile3)
% Convert chess coordinates to Cartesian
[x1, y1] = getCartesian(tile1);
[x2, y2] = getCartesian(tile2);
[x3, y3] = getCartesian(tile3);

% Determine direction of movement
direction1 = [x2-x1, y2-y1];
direction2 = [x3-x2, y3-y2];

% A turn occurs if the direction changes
isTurn = ~isequal(direction1, direction2);
end

function [xCoord, yCoord] = getCartesian(coord)
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
xCoord = find(alphabet == coord(1)); % MATLAB index starts at 1, so subtract 1
yCoord = str2double(coord(2)); % Convert second character to double for y-coordinate
end

function metrics_mz = ProcessAndAppendPathErrors(metrics_mz, mazeType, errorTypesInTrial, pathType, trialInBlock)
ignoreErrors = false;
tempErrorCounts = struct('ruleBreakingErrors', 0, ...
    'ruleAbidingErrors', 0, ...
    'backTrackErrors', 0, ...
    'perseverativeErrors', 0, ...
    'retouchErrors', 0);

for iError = 1:length(errorTypesInTrial)
    errorType = errorTypesInTrial{iError};

    % Reset error ignoring upon retouchCorrect
    if strcmp(errorType, 'retouchCorrect')
        ignoreErrors = false;
        continue;
    end

    % Skip counting errors until retouchCorrect is encountered
    if ignoreErrors
        continue;
    end

    % Count initial error and specific types
    if any(strcmp(errorType, {'ruleAbidingError', 'ruleBreakingError', 'backTrackError'}))
        tempErrorCounts.([errorType 's']) = tempErrorCounts.([errorType 's']) + 1; % Note the 's' to match struct field names
        ignoreErrors = true; % Start ignoring subsequent errors until retouchCorrect
    elseif contains(errorType, 'perseverative')
        tempErrorCounts.perseverativeErrors = tempErrorCounts.perseverativeErrors + 1;
        ignoreErrors = true; % Similar treatment for perseverative errors
    elseif strcmp(errorType, 'retouchError')
        tempErrorCounts.retouchErrors = tempErrorCounts.retouchErrors + 1;
        % Note: retouchError does not trigger ignoreErrors
    end
end

% Now, append the processed error counts to metrics_mz using AppendSpatialErrors
metrics_mz = AppendSubtypeErrors(metrics_mz, mazeType, tempErrorCounts, pathType, trialInBlock);
end
function metrics_mz = AppendSubtypeErrors(metrics_mz, mazeType, subTypeFields, subType, blockNum)
fields = fieldnames(subTypeFields);
for i = 1:numel(fields)
    fieldName = fields{i};
    currentCount = subTypeFields.(fieldName);
    metrics_mz.(mazeType).(subType).(fieldName){blockNum} = [metrics_mz.(mazeType).(subType).(fieldName){blockNum}, currentCount];
end
end




function plot_AcrossSessionAnalysis_SortedByMazeType (metrics_mz, res)
% Define fields to analyze
fields = {'normalizedTotalErrors', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors'};

% Prepare for plotting on the same figure
figure;
hold on; % Allows multiple plots on the same figure

% Extract number of turns and lengths for x-axis labels
numConditions = length(metrics_mz.Aggregate.numTurns);
xLabels = cell(1, numConditions);
for i = 1:numConditions
    numTurns = metrics_mz.Aggregate.numTurns{i}(1); % Assuming the first element is representative
    numLength = metrics_mz.Aggregate.numLength{i}(1); % Assuming the first element is representative
    xLabels{i} = sprintf('T%d.L%d', numTurns, numLength);
end

% Colors or markers for each field for distinction
markers = {'o-', '+-', '*-', 'x-'};
colors = {'blue', 'red', 'green', 'black'};

for f = 1:length(fields)
    fieldName = fields{f};
    data = metrics_mz.Aggregate.(fieldName); % Extract the data for this metric

    % Initialize arrays to store the averages and SEMs for each trial count
    avgData = zeros(1, length(data));
    semData = zeros(1, length(data));

    % Calculate average and standard deviation for each condition
    for i = 1:length(data)
        if ~isempty(data{i})
            n = length(data{i});  % Sample size

            avgData(i) = mean(data{i});
            semData(i) = std(data{i}) / sqrt(n);
        end
    end

    % Plot the current metric with error bars
    errorbar(1:length(avgData), avgData, semData, markers{f}, 'Color', colors{f}, 'DisplayName', strrep(fieldName, 'normalized', 'Normalized '));
end

% Customize the plot
title([res(1).subjectName '- Average Normalized Errors Across Maze Condition']);
xlabel('Maze Condition');
ylabel('Average Normalized Error Rates (Normalized by Path Length in Trial)');
xlim([0.5, length(avgData) + 0.5]); % Adjust the x-axis limits to improve visibility
set(gca, 'XTick', 1:length(avgData), 'XTickLabel', xLabels, 'XTickLabelRotation', 45); % Use custom x-axis labels
legend('show'); % Show legend to distinguish between metrics
hold off; % No more plots on this figure
end

function plot_AcrossSessionAnalysis_SortedBySession(metrics_mz, res)
% Prepare for plotting on the same figure
figure;
hold on; % Allows multiple plots on the same figure

% Extract session dates for x-axis labels
sessionDates = metrics_mz.Maze1.sessionDate; % Assuming session dates are consistent across 'Maze1' and 'Maze1Repeat'
numSessions = length(sessionDates);
xLabels = cell(1, numSessions);
for i = 1:numSessions
    if ~isempty(sessionDates{i})
        xLabels{i} = datestr(sessionDates{i}, 'mm/dd'); % Format date as 'MM/DD'
    else
        xLabels{i} = ['Session ', num2str(i)]; % Fallback label if date is missing
    end
end

% Specific dates to draw vertical lines
% verticalLineDates = {'02/12', '02/19', '02/26', '03/04', '03/11', '03/18'};
verticalLineDates = {'04/01', '04/08'};

% Define fields to analyze
fields = {'normalizedTotalErrors', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors'};

% Colors or markers for each field for distinction
markers = {'o-', '+-', '*-', 'x-'};
colors = {'blue', 'red', 'green', 'black'};

for f = 1:length(fields)
    fieldName = fields{f};
    data = metrics_mz.Aggregate.(fieldName); % Extract the data for this metric

    % Initialize arrays to store the averages and SEMs for each trial count
    avgData = zeros(1, length(numSessions));
    semData = zeros(1, length(numSessions));

    % Calculate average and standard deviation for each session
    for i = 1:numSessions
        if ~isempty(data{i})
            n = length(data{i});  % Sample size
            avgData(i) = mean(data{i});
            semData(i) = std(data{i}) / sqrt(n);  % Calculate SEM
        end
    end

    % Plot the current metric with error bars
    errorbar(1:numSessions, avgData, semData, markers{f}, 'Color', colors{f}, 'DisplayName', strrep(fieldName, 'normalized', 'Normalized '));
end

% Draw vertical lines at specific dates
for i = 1:length(verticalLineDates)
    idx = find(strcmp(xLabels, verticalLineDates{i}));
    if ~isempty(idx)
        line([idx idx], ylim, 'Color', 'black', 'LineStyle', '--', 'HandleVisibility', 'off');
    end
end

% Customize the plot
title([res(1).subjectName ' - Average Normalized Errors Across Session Date']);
xlabel('Session Date');
ylabel('Average Normalized Error Rates (Normalized by Path Length in Trial)');
xlim([0.5, numSessions + 0.5]); % Adjust the x-axis limits to improve visibility
set(gca, 'XTick', 1:numSessions, 'XTickLabel', xLabels, 'XTickLabelRotation', 45); % Use custom x-axis labels for each session
legend('show'); % Show legend to distinguish between metrics
hold off; % No more plots on this figure

end

function plot_AcrossSessionAnalysis_SortedByWeekday(metrics_mz, res)
% Prepare for plotting on the same figure
figure;
hold on; % Allows multiple plots on the same figure

% Define days of the week to analyze
daysOfTheWeek = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
% Define fields to analyze
fields = {'normalizedTotalErrors', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors'};

% Colors or markers for each field for distinction
markers = {'o-', '+-', '*-', 'x-'};
colors = {'blue', 'red', 'green', 'black'};

for f = 1:length(fields)
    fieldName = fields{f};

    % Arrays to store the averages and SEMs for each data point that exists
    avgData = [];
    semData = [];
    dayIndices = []; % To keep track of which days have data

    % Calculate average and standard deviation for each weekday
    for iDay = 1:length(daysOfTheWeek)
        data = metrics_mz.Aggregate.(fieldName){iDay}; % Assuming data is stored in a cell array

        if ~isempty(data)
            n = length(data);  % Sample size
            avgData(end+1) = mean(data);
            semData(end+1) = std(data) / sqrt(n);  % Calculate SEM
            dayIndices(end+1) = iDay; % Keep track of the day index
        end
    end

    % Plot data for this field
    errorbar(dayIndices, avgData, semData, markers{f}, 'Color', colors{f}, 'DisplayName', strrep(fieldName, 'normalized', 'Normalized '));
end

% Customize the plot
title([res(1).subjectName '- Average Normalized Errors by Weekday']);
xlabel('Weekday');
ylabel('Average Normalized Error Rates (Normalized by Path Length in Trial)');
xlim([0.5, length(daysOfTheWeek) + 0.5]); % Adjust the x-axis limits to improve visibility
set(gca, 'XTick', 1:length(daysOfTheWeek), 'XTickLabel', daysOfTheWeek, 'XTickLabelRotation', 45); % Use custom x-axis labels for each weekday
legend('show'); % Show legend to distinguish between metrics
hold off; % No more plots on this figure
end

function plot_AcrossSessionAnalysis_LearningCurve(metrics_mz, res, uniqueTrialCounts)
% Prepare for plotting on the same figure
figure;
hold on; % Allows multiple plots on the same figure
title('Learning Curve for Aggregate Error Analysis'); % Title for the plot

% Define fields to analyze
fields = {'normalizedTotalErrors', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors', 'normalizedMazeDurations'};

% Colors or markers for each field for distinction
markers = {'o-', '+-', '*-', 'x-', '--'};
colors = {'blue', 'red', 'green', 'black', 'magenta'};

% Loop through each error field to plot
for f = 1:length(fields)
    fieldName = fields{f};

    % Initialize arrays to store the averages and SEMs for each trial count
    avgData = zeros(1, length(uniqueTrialCounts));
    semData = zeros(1, length(uniqueTrialCounts));

    % Calculate average and SEM for each trial count
    for iTrialInBlock = 1:length(uniqueTrialCounts)
        data = metrics_mz.Aggregate.(fieldName){iTrialInBlock}; % Assuming data is stored similarly

        if ~isempty(data)
            % Remove outliers based on the interquartile range (IQR)
            Q1 = quantile(data, 0.25);
            Q3 = quantile(data, 0.75);
            IQR = Q3 - Q1;
            lowerBound = Q1 - 1.5 * IQR;
            upperBound = Q3 + 1.5 * IQR;

            dataWithoutOutliers = data(data >= lowerBound & data <= upperBound);

            % Display removed data count
            disp(['Data removed for ' fieldName ': ' num2str(length(data) - length(dataWithoutOutliers))]);

            % Proceed if there are any data points left after removing outliers
            if ~isempty(dataWithoutOutliers)
                n = length(dataWithoutOutliers);  % Sample size
                avgData(iTrialInBlock) = mean(dataWithoutOutliers);
                semData(iTrialInBlock) = std(dataWithoutOutliers) / sqrt(n);  % Calculate SEM
            end
        end
    end

    % Plot the current metric with error bars using SEM
    errorbar(1:length(uniqueTrialCounts), avgData, semData, markers{f}, 'Color', colors{f}, 'DisplayName', strrep(fieldName, 'normalized', 'Normalized '));
end

% Configure x-axis to show only whole number ticks corresponding to trial counts
xticks(1:length(uniqueTrialCounts)); % Set x-axis ticks to cover all trial counts explicitly

% Add plot legends, labels, and grid
title([res(1).subjectName '- Average Normalized Errors and Maze Durations by Trial in Block (Outliers Removed)']);
legend show;
xlim([0.5, length(uniqueTrialCounts) + 0.5]);
xlabel('Trial Count In Block');
ylabel('Average Normalized Error Rates and Maze Durations (Normalized by Path Length in Trial)');
grid on;
hold off;
end

function plot_AcrossSessionAnalysis_TurnVsStraightErrorRate(metrics_mz, res)
% Define error types
errorTypes = {'ruleBreakingErrors', 'ruleAbidingErrors', 'perseverativeErrors', 'retouchErrors', 'backTrackErrors'};

% Initialize arrays to store averages and SEMs for each error type
avgTurnData = zeros(1, length(errorTypes));
semTurnData = zeros(1, length(errorTypes));
avgStraightData = zeros(1, length(errorTypes));
semStraightData = zeros(1, length(errorTypes));

% Loop through each error type
for f = 1:length(errorTypes)
    fieldName = errorTypes{f};

    % Collect data for turns and straights
    allTurnData = [metrics_mz.Aggregate.turn.(fieldName){:}];
    allStraightData = [metrics_mz.Aggregate.straight.(fieldName){:}];

    % Calculate average and SEM for turns and straights
    avgTurnData(f) = mean(allTurnData);
    semTurnData(f) = std(allTurnData) / sqrt(length(allTurnData));
    avgStraightData(f) = mean(allStraightData);
    semStraightData(f) = std(allStraightData) / sqrt(length(allStraightData));
end

% Create a bar graph with error bars
figure;
hold on;  % Ensure no previous plots are interfering

% Bar graph with grouped bars
barData = [avgTurnData; avgStraightData]';
hb = bar(barData, 'grouped');

% Add error bars
% Calculate the positions for error bars for each group
for i = 1:length(hb)
    x = hb(i).XEndPoints;  % Get the end points of each set of bars
    if i == 1
        errorbar(x, barData(:, i), semTurnData, 'k', 'linestyle', 'none', 'CapSize', 10);
    else
        errorbar(x, barData(:, i), semStraightData, 'k', 'linestyle', 'none', 'CapSize', 10);
    end
end

% Set labels and title for the plot
xticks(1:length(errorTypes));
xticklabels(errorTypes);
xlabel('Error Type');
ylabel('Average Error Rate');
title([res(1).subjectName ' - Average Error Rate (Turn vs. Straight)']);
legend({'Turn', 'Straight'}, 'Location', 'best');
hold off;
end

function plot_AcrossSessionAnalysis_TowardsVsAway(metrics_mz, res)
% Define error types
errorTypes = {'ruleBreakingErrors', 'ruleAbidingErrors', 'perseverativeErrors', 'retouchErrors', 'backTrackErrors'};

% Initialize arrays to store averages and SEMs for each error type
avgTowardsData = zeros(1, length(errorTypes));
semTowardsData = zeros(1, length(errorTypes));
avgAwayData = zeros(1, length(errorTypes));
semAwayData = zeros(1, length(errorTypes));

% Loop through each error type
for f = 1:length(errorTypes)
    fieldName = errorTypes{f};

    % Collect data for towards and away
    allTowardsData = [metrics_mz.Aggregate.towards.(fieldName){:}];
    allAwayData = [metrics_mz.Aggregate.away.(fieldName){:}];

    % Calculate average and SEM for towards and away
    avgTowardsData(f) = mean(allTowardsData);
    semTowardsData(f) = std(allTowardsData) / sqrt(length(allTowardsData));
    avgAwayData(f) = mean(allAwayData);
    semAwayData(f) = std(allAwayData) / sqrt(length(allAwayData));
end

% Create a bar graph with error bars
figure;
hold on;  % Ensure no previous plots are interfering

% Bar graph with grouped bars
barData = [avgTowardsData; avgAwayData]';
hb = bar(barData, 'grouped');

% Add error bars
% Calculate the positions for error bars for each group
for i = 1:length(hb)
    x = hb(i).XEndPoints;  % Get the end points of each set of bars
    if i == 1
        errorbar(x, barData(:, i), semTowardsData, 'k', 'linestyle', 'none', 'CapSize', 10);
    else
        errorbar(x, barData(:, i), semAwayData, 'k', 'linestyle', 'none', 'CapSize', 10);
    end
end

% Set labels and title for the plot
xticks(1:length(errorTypes));
xticklabels(errorTypes);
xlabel('Error Type');
ylabel('Average Error Rate');
title([res(1).subjectName ' - Average Error Rate (Towards vs. Away)']);
legend({'Towards', 'Away'}, 'Location', 'best');
hold off;
end

function plot_AcrossSessionAnalysis_EarlyVsLate(metrics_mz, res)
% Define mazes and fields to analyze
mazes = {'Maze1', 'Maze1Repeat', 'Maze2'};
fields = {'normalizedTotalErrors', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors'};

% Colors for Early, Repeat, Late
colors = {[0.8500, 0.3250, 0.0980], [0, 0.4470, 0.7410], [0.4660, 0.6740, 0.1880]}; % MATLAB default color order
stageLabels = {'Early', 'Repeat', 'Late', }; % This will be used for the legend

figure;
hold on;

nGroups = length(fields);
nBars = length(mazes); % Adjusted for stages
groupWidth = 0.8; % Adjust if needed for wider group spacing
interGroupSpacing = 1.5;

for iField = 1:nGroups
    fieldName = fields{iField};

    for iMaze = 1:length(mazes)
        mazeName = mazes{iMaze};
        data = [metrics_mz.(mazeName).(fieldName){:}];

        % Calculate mean and SEM if stageData is not empty
        if ~isempty(data)
            mazeAverages(iMaze, iField) = mean(data);
            mazeSEMs(iMaze, iField) = std(data) / sqrt(length(data));
        end

        % Adjust xPosition calculation for better spacing between bars within each group
        xPosition = (iField - 1) * interGroupSpacing + (iMaze - 1) * (groupWidth / nBars) + 1;

        % Plotting the bar
        bar(xPosition, mazeAverages(iMaze, iField), groupWidth / nBars, 'FaceColor', colors{iMaze});

        % Overlaying error bars
        errorbar(xPosition, mazeAverages(iMaze, iField), mazeSEMs(iMaze, iField), 'k', 'linestyle', 'none');
    end
end

% Recalculate XTick positions to align with the center of each group
XTickPositions = 1:interGroupSpacing:nGroups*interGroupSpacing;

% Adjusting X-ticks and labels for clearer titles
set(gca, 'XTick', XTickPositions, 'XTickLabel', fields, 'XTickLabelRotation', 30);

% Create dummy bars for the legend for each stage only
for i = 1:length(colors)
    b(i) = bar(nan, 'FaceColor', colors{i}); % Dummy bar for legend
end
legend(b, stageLabels, 'Location', 'Best');
title([res(1).subjectName ' - Average Normalized Errors For Early, Late, and Repeat Mazes']);
ylabel('Average Normalized Error Rates (Normalized by Path Length in Trial)');

hold off;
end
function plot_AcrossSessionAnalysis_ReactionTimes(metrics_mz, res)
    % Define categories of reaction times to analyze
    categories = {'correctToCorrect', 'correctToError', 'errorToRetouch', ...
                  'retouchToCorrect', 'retouchToError', 'errorToError'};

    % Initialize arrays to store averages and SEMs for each category
    avgData = zeros(1, length(categories));
    semData = zeros(1, length(categories));

    % Define a color matrix with unique colors for each category
    colors = lines(length(categories)); % MATLAB's 'lines' colormap provides distinct colors

    % Loop through each category to collect and calculate data for Aggregate only
    for j = 1:length(categories)
        % Collect all reaction times for the current category under Aggregate
        reactionTimes = [metrics_mz.Aggregate.(categories{j}){:}];
        
        % Calculate average and SEM if there are any data points
        if ~isempty(reactionTimes)
            avgData(j) = mean(reactionTimes);
            semData(j) = std(reactionTimes) / sqrt(length(reactionTimes));
        end
    end

    % Create a bar graph with error bars
    figure;
    hold on;  % Ensure no previous plots are interfering

    % Bar graph for aggregate data with different colors for each bar
    for j = 1:length(categories)
        hb = bar(j, avgData(j), 'FaceColor', colors(j, :)); % Assign color to each bar
    end

    ylabel('Average Reaction Time (s)');
    title([res(1).subjectName ' - Aggregate Average Reaction Times by Selection Classifications']);
    set(gca, 'XTick', 1:length(categories), 'XTickLabel', categories);
    xtickangle(30);  % Adjust tick angle for better visibility if needed

    % Add error bars correctly positioned over each bar
    errorbar(1:length(categories), avgData, semData, 'k', 'linestyle', 'none', 'CapSize', 10);

    hold off;
end



function currentSessionData = AppendCurrentSessionData(currentTrialData, currentSessionData, blockDef)
iT = 0;
for jj = 1:length(currentTrialData)
    if currentTrialData(jj).AbortCode ~= 0
        continue; % Skip to the next iteration of the loop
    end
    iT = iT + 1;
    mazeStruct = jsondecode(blockDef(currentTrialData(jj).BlockCount).MazeDef);

    nTurns = mazeStruct.mNumTurns;
    mLength = mazeStruct.mNumSquares;

    currentSessionData.subjectName{iT} = currentTrialData(jj).SubjectID;%  trialData(iT).SubjectNum;
    currentSessionData.dayOfTheWeek(iT) = day(currentSessionData.sessionDatesFormatted, "longname");
    currentSessionData.mazeDefName{iT} = mazeStruct.mName;
    currentSessionData.mazePath{iT} = mazeStruct.mPath;


    currentSessionData.mazeTurnsLength(iT,1:2) = [nTurns, mLength];
    currentSessionData.blockNum(iT) = currentTrialData(jj).BlockCount;
    currentSessionData.blockName(iT) = {blockDef(currentTrialData(jj).BlockCount).BlockName};

    currentSessionData.trialInBlock(iT) = currentTrialData(jj).TrialCount_InBlock;
    currentSessionData.mazeDuration(iT) = currentTrialData(jj).MazeDuration;
    currentSessionData.sliderBarFilled(iT) = strcmp(currentTrialData(jj).SliderBarFilled,'True');
    currentSessionData.totalErrors(iT) = currentTrialData(jj).TotalErrors;
    currentSessionData.selectedTiles{iT} = strsplit(currentTrialData(jj).SelectedTiles, ',');
    [currentSessionData.errorTypes{iT}, currentSessionData.tileSpecificErrorCounts{iT}]  = determineSelectionError( currentSessionData.selectedTiles{iT}, mazeStruct.mPath);
    if isfield(currentTrialData(jj), 'ReactionTimePerSelectedTiles')
        currentSessionData.reactionTimes{iT} = strsplit(currentTrialData(jj).ReactionTimePerSelectedTiles, ',');
    end
    currentSessionData.correctTouches(iT) = currentTrialData(jj).CorrectTouches;
    currentSessionData.retouchCorrect(iT) = currentTrialData(jj).RetouchCorrect;
    currentSessionData.retouchErroneous(iT) = currentTrialData(jj).RetouchErroneous;
    currentSessionData.backTrackingErrors(iT) = currentTrialData(jj).BacktrackingErrors;
    currentSessionData.ruleAbidingErrors(iT) = currentTrialData(jj).Rule_AbidingErrors;
    currentSessionData.ruleBreakingErrors(iT) = currentTrialData(jj).Rule_BreakingErrors;
    currentSessionData.perseverativeRetouchErroneous(iT) = currentTrialData(jj).PerseverativeRetouchErrors;
    currentSessionData.perseverativeBackTrackingErrors(iT) = currentTrialData(jj).PerseverativeBackTrackErrors;
    currentSessionData.perseverativeRuleAbidingErrors(iT) = currentTrialData(jj).PerseverativeRuleAbidingErrors;
    currentSessionData.perseverativeRuleBreakingErrors(iT) = currentTrialData(jj).PerseverativeRuleBreakingErrors;

end
end
function metrics_mz = ProcessSpatialMazeData(metrics_mz, currentMazedata, mazeType, do_AcrossSessionAnalysis_TurnVsStraightErrorRate, do_AcrossSessionAnalysis_TowardsVsAwayEndTile)
for iTrial = 1:length(currentMazedata.trialInBlock)
    mPath = currentMazedata.mazePath{iTrial};
    turns = find(getTurnsAlongPath(mPath)) + 1;
    errorTypesInTrial = currentMazedata.errorTypes{iTrial};
    pathProgress = find(strcmp(errorTypesInTrial, 'correct'));
    pathType = '';
    directionType = '';

    if length(pathProgress) < 2
        continue;  % Need at least 1 'correct' to define towards/away and need atleast 2 'correct' to define turn/straight
    end

    startIndex = pathProgress(1) + 1;  % Initial error checking starts after the first correct
    for iTile = 2:length(pathProgress) % iTile represents the next tile along the path
        if iTile > length(mPath) || pathProgress(iTile) > length(errorTypesInTrial)
            break;  % Avoid indexing beyond the bounds
        end

        endIndex = pathProgress(iTile) - 1;  % Errors to be considered end before the next 'correct'
        errorTypesDuringSelection = errorTypesInTrial(startIndex:endIndex);

        currentTilePos = mPath{iTile - 1};
        targetTilePos = mPath{iTile};
        endTilePos = mPath{end};
        if iTile > 2
            if any(turns == iTile)
                pathType = 'turn';
            else
                pathType = 'straight';
            end
        end

        if(isMovingTowardsEnd(currentTilePos, targetTilePos, endTilePos))
            directionType = 'towards';
        else
            directionType = 'away';
        end

        blockNum = currentMazedata.blockNum(iTrial);
        startIndex = pathProgress(iTile) + 1;  % Update startIndex for the next set of errors

        if(~strcmp(pathType, '') && (do_AcrossSessionAnalysis_TurnVsStraightErrorRate == 1))
            metrics_mz = ProcessAndAppendPathErrors(metrics_mz, mazeType, errorTypesDuringSelection, pathType, blockNum);
            pathType = '';
        end
        if(~strcmp(directionType, '') && (do_AcrossSessionAnalysis_TowardsVsAwayEndTile == 1))
            metrics_mz = ProcessAndAppendPathErrors(metrics_mz, mazeType, errorTypesDuringSelection, directionType, blockNum);
            directionType = '';
        end

    end
end
end

function metrics_mz = ProcessReactionTimeData(metrics_mz, currentMazeData, mazeType)

for iTrial = 1 : length(currentMazeData.trialInBlock)
    errorTypes = currentMazeData.errorTypes{iTrial};
    reactionTimes = currentMazeData.reactionTimes{iTrial};

    
    % Collect reaction times based on error transitions
    for iError = 2:length(errorTypes)-1
        currentType = errorTypes{iError -1};
        nextType = errorTypes{iError};
currentReactionTime = str2double(reactionTimes{iError}); % Convert to double if necessary

            % Check for NaN in case the conversion fails
            if isnan(currentReactionTime)
                continue; % Skip this iteration if the reaction time isn't a valid number
            end
        if strcmp(currentType, 'correct') && strcmp(nextType, 'correct')
            metrics_mz.(mazeType).correctToCorrect = [metrics_mz.(mazeType).correctToCorrect, currentReactionTime];
        elseif strcmp(currentType, 'correct') && contains(lower(nextType), 'error')
            metrics_mz.(mazeType).correctToError = [metrics_mz.(mazeType).correctToError, currentReactionTime];
        elseif contains(lower(currentType), 'error') && contains(lower(nextType), 'error')
            metrics_mz.(mazeType).errorToError = [metrics_mz.(mazeType).correctToError, currentReactionTime];
        elseif contains(lower(currentType), 'error') && strcmp(nextType, 'retouchCorrect')
            metrics_mz.(mazeType).errorToRetouch = [metrics_mz.(mazeType).errorToRetouch, currentReactionTime];
        elseif strcmp(currentType, 'retouchCorrect') && strcmp(nextType, 'correct')
            metrics_mz.(mazeType).retouchToCorrect = [metrics_mz.(mazeType).retouchToCorrect, currentReactionTime];
        elseif strcmp(currentType, 'retouchCorrect') && contains(lower(nextType), 'error')
            metrics_mz.(mazeType).retouchToError = [metrics_mz.(mazeType).retouchToError, currentReactionTime];
        end
    end
end

end

