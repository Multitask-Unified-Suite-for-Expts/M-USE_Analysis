
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
do_AcrossSessionAnalysis_LearningCurve_ErrorRateDifferences = 0;


do_WithinSessionAnalysis = 0; % NEEDS TO BE UPDATED

DO_SAVEFIGURES = 0;

PLOT_INDIVIDUALSESSIONS = 0;

MIN_TURNS = 4; % SET TO 0 TO IGNORE TRIAL DATA FILTERING
MIN_LENGTH = 14;
MAX_TURNS = 0;
MAX_LENGTH = 0;

% --- --- --- --- --- --- --- ---
% --- folder with preprocecced data
% --- --- --- --- --- --- --- ---
HOME_FOLDER = ['C:\Users\Sorti\OneDrive\Documents\GitHub\M-USE_Analysis\M-USE DataProcessingPipeline' filesep];
%'/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS'
MUSEMATFOLDERNames = {}; iResultFolder = ''; iResultFile = '';

 SessionID_MZ = {'Frey_MZ_all_01'};
% SessionID_MZ = {'Wotan_MZ_all_01'};

% MUSEMATFOLDERNames = 'MUSEMAT01_VS_EC_FL_MZG_Frey'
switch SessionID_MZ{1}
    case {'Frey_MZ_all_01'}
        MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_FL_MZG_Frey'
        % MUSEMATFOLDERNames{2} = 'MUSEMAT01_VS_EC_FL_MRT_Frey'
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Frey_01_20240322'; %'WM_01_20231205';%'WM_01_20230907';
    case {'Wotan_MZ_all_01'}
        MUSEMATFOLDERNames{1} =  'MUSEMAT01_WM_EC_FL_MZG_Wotan';
        %MUSEMATFOLDERNames{2} =  'MUSEMAT01_VS_EC_FL_MRT_Wotan';
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Wotan_01_20240322';%'WM_01_20230907';
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

                    trialData = in.dat.trialData{idx};
                    blockData = in.dat.blockData{idx};
                    blockDef = table2struct(in.dat.cfg_BlockDef{idx});
                    res(iD).subjectName = trialData(1).SubjectID;
                    iT = 0;

                    sessionDates = cellfun(@(x) regexp(x, 'Session_(\d+_\d+)', 'tokens'), {res(iD).dataset}, 'UniformOutput', false);
                    sessionDates = cellfun(@(x) x{1}{1}, sessionDates, 'UniformOutput', false); % Convert tokens to strings

                    currentSessionData = struct(); % data for the current session
                    currentSessionData.sessionDatesFormatted = datetime(sessionDates, 'InputFormat', 'MM_dd', 'Format', 'MM/dd'); % Convert to datetime
                    res(iD).sessionDatesFormatted = currentSessionData.sessionDatesFormatted;

                    for jj = 1:length(trialData)
                        if trialData(jj).AbortCode ~= 0
                            continue; % Skip to the next iteration of the loop
                        end

                        mazeStruct = jsondecode(blockDef(trialData(jj).BlockCount).MazeDef);

                        nTurns = mazeStruct.mNumTurns;
                        mLength = mazeStruct.mNumSquares;

                        iT = iT + 1;
                        currentSessionData.subjectName{iT} = trialData(jj).SubjectID;%  trialData(iT).SubjectNum;
                        currentSessionData.dayOfTheWeek(iT) = day(currentSessionData.sessionDatesFormatted, "longname");
                        currentSessionData.mazeDefName{iT} = mazeStruct.mName;


                        currentSessionData.mazeTurnsLength(iT,1:2) = [nTurns, mLength];
                        currentSessionData.blockNum(iT) = trialData(jj).BlockCount;
                        currentSessionData.trialInBlock(iT) = trialData(jj).TrialCount_InBlock;
                        currentSessionData.mazeDuration(iT) = trialData(jj).MazeDuration;
                        currentSessionData.sliderBarFilled(iT) = strcmp(trialData(jj).SliderBarFilled,'True');
                        currentSessionData.totalErrors(iT) = trialData(jj).TotalErrors;
                        currentSessionData.selectedTiles{iT} = strsplit(trialData(jj).SelectedTiles, ',');
                        currentSessionData.correctTouches(iT) = trialData(jj).CorrectTouches;
                        currentSessionData.retouchCorrect(iT) = trialData(jj).RetouchCorrect;
                        currentSessionData.retouchErroneous(iT) = trialData(jj).RetouchErroneous;
                        currentSessionData.backTrackingErrors(iT) = trialData(jj).BacktrackingErrors;
                        currentSessionData.ruleAbidingErrors(iT) = trialData(jj).Rule_AbidingErrors;
                        currentSessionData.ruleBreakingErrors(iT) = trialData(jj).Rule_BreakingErrors;
                        currentSessionData.perseverativeRetouchErroneous(iT) = trialData(jj).PerseverativeRetouchErrors;
                        currentSessionData.perseverativeBackTrackingErrors(iT) = trialData(jj).PerseverativeBackTrackErrors;
                        currentSessionData.perseverativeRuleAbidingErrors(iT) = trialData(jj).PerseverativeRuleAbidingErrors;
                        currentSessionData.perseverativeRuleBreakingErrors(iT) = trialData(jj).PerseverativeRuleBreakingErrors;


                        % Save the collected data back to the main structure
                        res(iD).data.(strrep(label{1}, '_', '')) = currentSessionData;
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
        end

        uniqueMazeConfigurations = unique(allMazes, 'rows');
    end

    if (do_AcrossSessionAnalysis_LearningCurve == 1 || do_AcrossSessionAnalysis_LearningCurve_ErrorRateDifferences == 1 || do_WithinSessionAnalysis == 1)
        allTrials= [];
        for i = 1:length(res)
            if isfield(res(i).data, 'Maze1') && ~isempty(fieldnames(res(i).data.Maze1))
                allTrials = [allTrials; res(i).data.Maze1.trialInBlock'];
            end
            if isfield(res(i).data, 'Maze1Repeat') && ~isempty(fieldnames(res(i).data.Maze1Repeat))
                allTrials = [allTrials; res(i).data.Maze1Repeat.trialInBlock'];
            end
        end

        uniqueTrialCounts = unique(allTrials, 'rows');
    end


    % Initialize the empty arrays to store data for the respective analysis
    % type
    mazeTypes = {'Maze1', 'Maze1Repeat', 'Aggregate'};
    metrics = {'normalizedTotalErrors', 'normalizedDurations', ...
        'normalizedRuleBreakingErrors', 'normalizedRuleAbidingErrors', ...
        'normalizedPerseverativeErrors'};

    if do_AcrossSessionAnalysis_SortedByMazeType == 1
        metrics = [metrics, {'numTurns', 'numLength'}];
        metricsLength = size(uniqueMazeConfigurations, 1);
    elseif do_AcrossSessionAnalysis_SortedBySession == 1
        metrics = [metrics, {'sessionDate'}];
        metricsLength = length(res);
    elseif do_AcrossSessionAnalysis_SortedByWeekday == 1
        metrics = [metrics, {'dayOfTheWeek'}];
        metricsLength = 5;
    elseif (do_AcrossSessionAnalysis_LearningCurve == 1 || do_WithinSessionAnalysis == 1)
        metrics = [metrics, {'trialInBlock'}];
        metricsLength = size(uniqueTrialCounts, 1);
    elseif(do_AcrossSessionAnalysis_LearningCurve_ErrorRateDifferences == 1 )
        metrics = [metrics, {'trialInBlock', 'errorRateDifference'}];
        metricsLength = size(uniqueTrialCounts, 1);
    end


    % Initialize the structure
    metrics_mz = struct();

    % Iterate through each maze type and metric to initialize arrays

    for i = 1:length(mazeTypes)
        mazeType = mazeTypes{i};
        for j = 1:length(metrics)
            metric = metrics{j};
            metrics_mz.(mazeType).(metric) = cell(metricsLength, 1);
        end
    end

    if(do_WithinSessionAnalysis ==1 )
        SESSION_DATE = "03/18";
        sessionIdx = find(res.sessionDatesFormatted == SESSION_DATE);

        initialMazeData = res(sessionIdx).data.Maze1;
        repeatMazeData = res(sessionIdx).data.Maze1Repeat;
        metrics_mz = AppendMetrics(metrics_mz, 'Maze1', initialMazeData, iL, filteredMaze1Indices);
        metrics_mz.Maze1.sessionDate{sessionIdx} = initialMazeData.sessionDatesFormatted;

        metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMazeData, iL, filteredMaze1Indices);

        metrics_mz.Maze1Repeat.sessionDate{sessionIdx} = repeatMazeData.sessionDatesFormatted;

    else
        for sessionIdx = 1:length(res)
            initialMazeData = res(sessionIdx).data.Maze1;
            repeatMazeData = res(sessionIdx).data.Maze1Repeat;
            if do_AcrossSessionAnalysis_SortedByMazeType == 1
                for iL = 1:size(uniqueMazeConfigurations, 1)
                    filteredMaze1Indices = find(initialMazeData.mazeTurnsLength(:,1) ==  uniqueMazeConfigurations(iL,1) & initialMazeData.mazeTurnsLength(:,2) == uniqueMazeConfigurations(iL,2));
                    % if ~isempty(filteredMaze1Indices)
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1', initialMazeData, iL, filteredMaze1Indices);
                    metrics_mz.Maze1.numTurns{iL} = uniqueMazeConfigurations(iL,1);
                    metrics_mz.Maze1.numLength{iL} = uniqueMazeConfigurations(iL,2);

                    filteredMaze1RepeatIndices = find(repeatMazeData.mazeTurnsLength(:,1) ==  uniqueMazeConfigurations(iL,1) & repeatMazeData.mazeTurnsLength(:,2) == uniqueMazeConfigurations(iL,2));
                    % if ~isempty(filteredMaze1RepeatIndices)
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMazeData, iL, filteredMaze1RepeatIndices);
                    metrics_mz.Maze1Repeat.numTurns{iL} = uniqueMazeConfigurations(iL,1);
                    metrics_mz.Maze1Repeat.numLength{iL} = uniqueMazeConfigurations(iL,2);
                end

            elseif do_AcrossSessionAnalysis_SortedBySession == 1
                metrics_mz = AppendMetrics(metrics_mz, 'Maze1', initialMazeData, sessionIdx, 1:height(initialMazeData.mazeTurnsLength));
                metrics_mz.Maze1.sessionDate{sessionIdx} = initialMazeData.sessionDatesFormatted;

                metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMazeData, sessionIdx, 1:height(repeatMazeData.mazeTurnsLength));
                metrics_mz.Maze1Repeat.sessionDate{sessionIdx} = repeatMazeData.sessionDatesFormatted;

            elseif do_AcrossSessionAnalysis_SortedByWeekday == 1
                daysOfTheWeek = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
                for iDay = 1:5
                    day = daysOfTheWeek(iDay);
                    filteredMaze1Indices = find(strcmpi(initialMazeData.dayOfTheWeek(1,:), day));
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1', initialMazeData, iDay, filteredMaze1Indices);
                    metrics_mz.Maze1.dayOfTheWeek(iDay) = day;

                    filteredMaze1RepeatIndices = find(strcmpi(repeatMazeData.dayOfTheWeek(1,:), day));
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMazeData, iDay, filteredMaze1RepeatIndices);
                    metrics_mz.Maze1Repeat.dayOfTheWeek(iDay) = day;
                end
            elseif (do_AcrossSessionAnalysis_LearningCurve == 1)
                for iTrialInBlock = 1:height(uniqueTrialCounts)
                    filteredMaze1Indices = find(initialMazeData.trialInBlock == iTrialInBlock);
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1', initialMazeData, iTrialInBlock, filteredMaze1Indices);
                    metrics_mz.Maze1.trialInBlock = iTrialInBlock;

                    filteredMaze1RepeatIndices = find(repeatMazeData.trialInBlock == iTrialInBlock);
                    metrics_mz = AppendMetrics(metrics_mz, 'Maze1Repeat', repeatMazeData, iTrialInBlock, filteredMaze1RepeatIndices);
                    metrics_mz.Maze1Repeat.trialInBlock = iTrialInBlock;
                end

            elseif (do_AcrossSessionAnalysis_LearningCurve_ErrorRateDifferences == 1)
                % Process for initialMazeData and repeatMazeData
                mazeDataArray = {initialMazeData; repeatMazeData};
                for i = 1:length(mazeDataArray)
                    currentMazeData = mazeDataArray{i};
                    uniqueBlocks = unique(currentMazeData.blockNum);

                    for iBlock = 1:length(uniqueBlocks)
                        blockNum = uniqueBlocks(iBlock);

                        % Find all trials in this block, sorted
                        blockTrialsIndices = find(currentMazeData.blockNum == blockNum);
                        [~, sortIdx] = sort(currentMazeData.trialInBlock(blockTrialsIndices));
                        sortedTrialsIndices = blockTrialsIndices(sortIdx);

                        prevError = NaN; % Initialize as NaN for the first trial in each block

                        for iTrialIdx = 1:length(sortedTrialsIndices)
                            idx = sortedTrialsIndices(iTrialIdx);
                            trialInBlock = currentMazeData.trialInBlock(idx);
                            proportionalError = currentMazeData.normalizedError(idx); % Example calculation

                            % Calculate errorRateDifference
                            if ~isnan(prevError) % Skip for the first trial in the block
                                errorRateDifference = proportionalError - prevError;
                                % Append this errorRateDifference to metrics_mz
                                metrics_mz.(mazeType{1}).errorRateDifference{trialInBlock} = [metrics_mz.(mazeType{1}).errorRateDifference{trialInBlock}; errorRateDifference];
                            else
                                % For the first trial in the block, append NaN
                                metrics_mz.(mazeType{1}).errorRateDifference{trialInBlock} = [metrics_mz.(mazeType{1}).errorRateDifference{trialInBlock}; NaN];
                            end

                            % Append proportionalError for the current trial
                            metrics_mz.(mazeType{1}).normalizedTotalErrors{trialInBlock} = [metrics_mz.(mazeType{1}).normalizedTotalErrors{trialInBlock}; proportionalError];

                            % Update prevError for the next iteration
                            prevError = proportionalError;
                        end
                    end
                end
            end

        end

        metrics_mz = AggregateMetrics(metrics_mz);

    end
end

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

        figure('Color', 'w');
        set(gcf, 'Position', [100, 100, 1200, 600]); % Adjust figure size to fit 4 plots in a row
        nc = 2; % Number of columns in subplot grid
        nr = 1; % Number of rows in subplot grid

        % Preparing labels for x-axis
        mazeConditions = strings(1, height(metrics_mz.numTurnsTiles));
        for j = 1:height(metrics_mz.numTurnsTiles)
            mazeConditions(j) = "T" + num2str(metrics_mz.numTurnsTiles(j,1)) + "/L" + num2str(metrics_mz.numTurnsTiles(j,2));
        end

        % Function to calculate SEM
        calcSEM = @(data) std(data, 'omitnan') / sqrt(length(data));

        % Iterate over selected mazes for plotting
        subplot(nr, nc, 1);

        hold on;

        % Initialize arrays for legend handles and their labels
        h = zeros(3, 1);
        legendLabels = {'First Run', 'Repeat Run', 'Aggregate'};

        for j = 1:height(metrics_mz.numTurnsTiles)
            selectedMazeIndex = j;
            % Aggregated errors for comparison
            durationsFirst = metrics_mz.normalizedDurations{selectedMazeIndex, 1};
            durationsRepeat = metrics_mz.normalizedDurations{selectedMazeIndex, 2};
            durationsAggregate = metrics_mz.normalizedDurations{selectedMazeIndex, 3};

            medianDFirst = median(durationsFirst, 'omitnan');
            semDFirst = calcSEM(durationsFirst);
            medianDRepeat = median(durationsRepeat, 'omitnan');
            semDRepeat = calcSEM(durationsRepeat);
            medianDAggregate = median(durationsAggregate, 'omitnan');
            semDAggregate = calcSEM(durationsAggregate);

            % Plot each with different color/marker for distinction
            errorbar(j-0.2, medianDFirst, semDFirst, 's', 'MarkerFaceColor', 'r', 'Color', 'r');
            errorbar(j, medianDRepeat, semDRepeat, 'd', 'MarkerFaceColor', 'b', 'Color', 'b');
            errorbar(j+0.2, medianDAggregate, semDAggregate, 'o', 'MarkerFaceColor', 'g', 'Color', 'g');

            % On the first iteration, save handles for the legend
            if j == 1
                h(1) = errorbar(NaN, NaN, 's', 'MarkerFaceColor', 'r', 'Color', 'r');
                h(2) = errorbar(NaN, NaN, 'd', 'MarkerFaceColor', 'b', 'Color', 'b');
                h(3) = errorbar(NaN, NaN, 'o', 'MarkerFaceColor', 'g', 'Color', 'g');
            end
        end

        xlim([0.5, height(metrics_mz.numTurnsTiles) + 0.5]);
        xticks(1:height(metrics_mz.numTurnsTiles));
        xticklabels(mazeConditions);
        xlabel('Maze Condition');
        ylabel('Normalized Duration');
        title('Normalized Durations - First vs. Repeat vs. Aggregate');

        % Use the handles for the legend
        legend(h, legendLabels, 'Location', 'best');

        hold off;

        % Comparing Total Errors - First vs. Repeat vs. Aggregate
        subplot(nr, nc, 2);
        hold on;

        % Initialize arrays for legend handles and their labels
        h = zeros(3, 1);
        legendLabels = {'First Run', 'Repeat Run', 'Aggregate'};

        for j = 1:height(metrics_mz.numTurnsTiles)
            selectedMazeIndex = j;
            % Aggregated errors for comparison
            totalErrorsFirst = metrics_mz.normalizedTotalErrors{selectedMazeIndex, 1};
            totalErrorsRepeat = metrics_mz.normalizedTotalErrors{selectedMazeIndex, 2};
            totalErrorsAggregate = metrics_mz.normalizedTotalErrors{selectedMazeIndex, 3};

            medianTEFirst = median(totalErrorsFirst, 'omitnan');
            semTEFirst = calcSEM(totalErrorsFirst);
            medianRepeat = median(totalErrorsRepeat, 'omitnan');
            semTERepeat = calcSEM(totalErrorsRepeat);
            medianTEAggregate = median(totalErrorsAggregate, 'omitnan');
            semTEAggregate = calcSEM(totalErrorsAggregate);

            % Plot each with different color/marker for distinction
            errorbar(j-0.2, medianTEFirst, semTEFirst, 's', 'MarkerFaceColor', 'r', 'Color', 'r');
            errorbar(j, medianRepeat, semTERepeat, 'd', 'MarkerFaceColor', 'b', 'Color', 'b');
            errorbar(j+0.2, medianTEAggregate, semTEAggregate, 'o', 'MarkerFaceColor', 'g', 'Color', 'g');

            % On the first iteration, save handles for the legend
            if j == 1
                h(1) = errorbar(NaN, NaN, 's', 'MarkerFaceColor', 'r', 'Color', 'r');
                h(2) = errorbar(NaN, NaN, 'd', 'MarkerFaceColor', 'b', 'Color', 'b');
                h(3) = errorbar(NaN, NaN, 'o', 'MarkerFaceColor', 'g', 'Color', 'g');
            end
        end

        xlim([0.5, height(metrics_mz.numTurnsTiles) + 0.5]);
        xticks(1:height(metrics_mz.numTurnsTiles));
        xticklabels(mazeConditions);
        xlabel('Maze Condition');
        title('Normalized Total Errors - First vs. Repeat vs. Aggregate');
        ylabel('Normalized Total Errors');

        % Use the handles for the legend
        legend(h, legendLabels, 'Location', 'best');

        hold off;


    end


    if do_AcrossSessionAnalysis_SortedBySession == 1
        % Extract session dates from dataset names

        % Initialize arrays to hold the means for plotting
        normalizedDurationsDiff = [];
        normalizedTotalErrorsDiff = [];
        normalizedRuleAbidingErrorsDiff = [];
        normalizedRuleBreakingErrorsDiff = [];
        normalizedPerseverationErrorsDiff = [];


        % Calculate the difference between the metrics for Maze1Repeat and Maze1
        for i = 1:length(metrics_mz1)
            normalizedDurationsDiff = [normalizedDurationsDiff; metrics_mz1Repeat(i).meanNormalizedDurations - metrics_mz1(i).meanNormalizedDurations];
            normalizedTotalErrorsDiff = [normalizedTotalErrorsDiff; metrics_mz1Repeat(i).meanNormalizedTotalErrors - metrics_mz1(i).meanNormalizedTotalErrors];
            normalizedRuleAbidingErrorsDiff = [normalizedRuleAbidingErrorsDiff; metrics_mz1Repeat(i).meanNormalizedRuleAbidingErrors - metrics_mz1(i).meanNormalizedRuleAbidingErrors];
            normalizedRuleBreakingErrorsDiff = [normalizedRuleBreakingErrorsDiff; metrics_mz1Repeat(i).meanNormalizedRuleBreakingErrors - metrics_mz1(i).meanNormalizedRuleBreakingErrors];
            normalizedPerseverationErrorsDiff = [normalizedPerseverationErrorsDiff; metrics_mz1Repeat(i).meanNormalizedPerseverationErrors - metrics_mz1(i).meanNormalizedPerseverationErrors];

        end

        % Create a figure for the plots
        figure;

        % Specify dates for vertical lines
        specificDates = ["02/12", "02/19", "02/26", "03/04", "03/11", "03/18"];
        year = 2024; % Assuming the year is 2024, adjust accordingly
        datesForLines = datetime(specificDates + "/" + num2str(year), 'InputFormat', 'MM/dd/yyyy');

        % Titles for each subplot
        plotTitles = {'Normalized Durations Difference', 'Normalized Total Errors Difference', 'Normalized Rule Abiding Errors Difference', 'Normalized Rule Breaking Errors Difference', 'Normalized Perseveration Errors Difference'};

        % Loop through each subplot to plot data
        for plotIndex = 1:5
            subplot(3, 2, plotIndex);
            hold on;

            % Select the appropriate dataset for each subplot
            switch plotIndex
                case 1
                    dataToPlot = normalizedDurationsDiff;
                case 2
                    dataToPlot = normalizedTotalErrorsDiff;
                case 3
                    dataToPlot = normalizedRuleAbidingErrorsDiff;
                case 4
                    dataToPlot = normalizedRuleBreakingErrorsDiff;
                case 5
                    dataToPlot = normalizedPerseverationErrorsDiff;
            end

            % Plot lines and scatter for the selected dataset
            plot(sessionDatesFormatted, dataToPlot, 'LineWidth', 1.5, 'DisplayName', 'Repeat - Initial');
            scatter(sessionDatesFormatted, dataToPlot, 'filled', 'HandleVisibility', 'off');

            % Draw vertical lines at specific dates
            for xValue = datesForLines
                xline(xValue, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
            end
            yline(0, 'r', 'LineWidth', 1, 'HandleVisibility', 'off');

            title(plotTitles{plotIndex});
            xlabel('Session Date');
            ylabel('Difference');
            xticks(sessionDatesFormatted);  % Set x-axis ticks at each date
            xtickformat('MM/dd');  % Format the tick labels
            xtickangle(60);  % Rotate labels for readability
            xlim([min(sessionDatesFormatted) max(sessionDatesFormatted)]);  % Adjust x-axis limits

            hold off;
        end

        % Adjust figure size and display
        set(gcf, 'Position', [50, 50, 1800, 1800]); % Adjust size as needed

        % Add legend in the last plot
        legend('Location', 'best');
    end
    if do_AcrossSessionAnalysis_SortedByWeekday == 1

        weekdays = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'};
        nDays = numel(weekdays); % Number of weekdays to plot
        metricsNames = {'normalizedDurations', 'normalizedTotalErrors', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors'};
        colors = {'b', 'r', 'g', 'm', 'c'}; % Assign a color for each metric for clarity

        % Function to calculate SEM
        calcSEM = @(data) std(data, 'omitnan') / sqrt(length(data));

        % Creating a figure to hold all subplots
        figure('Color', 'w', 'Position', [100, 100, 1200, 900]); % Adjust size as needed
        sgtitle([res(1).subjectName ' Weekday Performance'])

        for iMetric = 1:length(metricsNames)
            metricName = metricsNames{iMetric};
            color = colors{iMetric};

            % Create a subplot for each metric
            subplot(3, 2, iMetric); % Adjust the grid size if necessary
            hold on;

            for iDay = 1:nDays
                dayName = weekdays{iDay};

                % Assuming the existence of these metrics for each day
                metricData = metrics_mz.(dayName).(metricName);
                medianMetric = median(metricData, 'omitnan');
                semMetric = calcSEM(metricData);

                % Plot with error bars
                errorbar(iDay, medianMetric, semMetric, 'o', 'MarkerFaceColor', color, 'Color', color);
            end

            xlim([0.5, nDays + 0.5]);
            xticks(1:nDays);
            xticklabels(weekdays);
            xlabel('Day of the Week');
            ylabel(strrep(metricName, 'normalized', ''));
            title(strrep(metricName, 'normalized', 'Normalized '));

            hold off;
        end

    end
    if(do_AcrossSessionAnalysis_LearningCurve == 1)
        figure;

        % Improvement threshold
        improvementThreshold = 0.50;

        % --- Subplot 1: Scatter plot of raw normalizedTotalErrors ---
        subplot(2, 1, 1); % This specifies a 2-row, 1-column grid of subplots, and activates the first subplot.
        hold on;
        trialCounts = length(metrics_mz.Aggregate.normalizedTotalErrors);
        for iTrial = 1:trialCounts
            if ~isempty(metrics_mz.Aggregate.normalizedTotalErrors{iTrial})
                trialIndices = iTrial * ones(size(metrics_mz.Aggregate.normalizedTotalErrors{iTrial}));
                errorRates = metrics_mz.Aggregate.normalizedTotalErrors{iTrial};
                isBelowThreshold = errorRates <= improvementThreshold;
                scatter(trialIndices(~isBelowThreshold), errorRates(~isBelowThreshold), 'MarkerEdgeColor', 'b');
                scatter(trialIndices(isBelowThreshold), errorRates(isBelowThreshold), 'MarkerEdgeColor', 'r');
            end
        end
        xlabel('Trial Count in Block');
        ylabel('Normalized Total Error Rate');
        title('Raw Error Rates per Trial');
        legend('Above Threshold', 'Below Threshold', 'Location', 'best');
        grid on;
        hold off;

        % --- Subplot 2: Error bars for average normalizedTotalErrors per trialCountInBlock, excluding outliers ---
        subplot(2, 1, 2); % Activates the second subplot
        trialCounts = length(metrics_mz.Aggregate.normalizedTotalErrors);
        avgErrorRates = NaN(1, trialCounts);
        stderrErrorRates = NaN(1, trialCounts); % Standard Error of the mean, excluding outliers
        for iTrial = 1:trialCounts
            if ~isempty(metrics_mz.Aggregate.normalizedTotalErrors{iTrial})
                % Extract error rates for the current trial
                errorRates = metrics_mz.Aggregate.normalizedTotalErrors{iTrial};

                % Identify and remove outliers using the interquartile range (IQR)
                Q1 = quantile(errorRates, 0.25);
                Q3 = quantile(errorRates, 0.75);
                IQR = Q3 - Q1;
                lowerBound = Q1 - 1.5 * IQR;
                upperBound = Q3 + 1.5 * IQR;
                nonOutliers = errorRates(errorRates >= lowerBound & errorRates <= upperBound);
                disp(size(errorRates, 2) - length(nonOutliers))
                % Calculate the average and standard error for non-outlier error rates
                if ~isempty(nonOutliers) % Ensure there are non-outlier points to calculate
                    avgErrorRates(iTrial) = mean(nonOutliers);
                    stderrErrorRates(iTrial) = std(nonOutliers) / sqrt(length(nonOutliers));
                end
            end
        end

        % Plot the error bars for average error rates excluding outliers
        errorbar(1:trialCounts, avgErrorRates, stderrErrorRates, 'o-');
        xlabel('Trial Count in Block');
        ylabel('Average Normalized Total Error Rate ');
        title('Average Normalized Total Error Rate per Trial with Standard Error (Excluding Outliers)');
        grid on;




        % metricsNames = {'normalizedTotalErrors', 'normalizedDurations', 'normalizedRuleAbidingErrors', 'normalizedRuleBreakingErrors', 'normalizedPerseverativeErrors'};
        %
        % % Colors assigned for Maze1, Maze1Repeat, and Aggregate
        % mazeColors = {'b', 'r', 'k'}; % Blue for Maze1, Red for Maze1Repeat, Black for Aggregate
        %
        % % Function to calculate SEM
        % calcSEM = @(data) std(data, 'omitnan') / sqrt(length(data));
        %
        % % Creating a figure to hold all subplots
        % figure('Color', 'w', 'Position', [100, 100, 1200, 900]); % Adjust size as needed
        % % Construct the title string
        % titleStr = sprintf('%s Learning Curve', res(1).subjectName);
        %
        % % Set the title using sgtitle without any LaTeX formatting issues
        % sgtitle(titleStr, 'Interpreter', 'none');
        % uniqueTrials = unique([metrics_mz.Aggregate.trialInBlock]); % Assume Aggregate covers the same trials
        % nTrials = numel(uniqueTrials); % Number of unique trial blocks to plot
        %
        % for iMetric = 1:length(metricsNames)
        %     metricName = metricsNames{iMetric};
        %
        %     % Create a subplot for each metric
        %     subplot(ceil(length(metricsNames) / 2), 2, iMetric); % Adjust grid size if necessary
        %     hold on;
        %
        %     % % Plotting for Maze1
        %     % means = arrayfun(@(t) mean(metrics_mz.Maze1.(metricName)(metrics_mz.Maze1.trialInBlock == t), 'omitnan'), uniqueTrials);
        %     % sems = arrayfun(@(t) calcSEM(metrics_mz.Maze1.(metricName)(metrics_mz.Maze1.trialInBlock == t)), uniqueTrials);
        %     % errorbar(uniqueTrials, means, sems, 'o-', 'Color', mazeColors{1}, 'DisplayName', 'Maze1', 'MarkerFaceColor', mazeColors{1});
        %     %
        %     % % Plotting for Maze1Repeat
        %     % means = arrayfun(@(t) mean(metrics_mz.Maze1Repeat.(metricName)(metrics_mz.Maze1Repeat.trialInBlock == t), 'omitnan'), uniqueTrials);
        %     % sems = arrayfun(@(t) calcSEM(metrics_mz.Maze1Repeat.(metricName)(metrics_mz.Maze1Repeat.trialInBlock == t)), uniqueTrials);
        %     % errorbar(uniqueTrials, means, sems, 'o-', 'Color', mazeColors{2}, 'DisplayName', 'Maze1Repeat', 'MarkerFaceColor', mazeColors{2});
        %
        %     % Plotting for Aggregate
        %     means = arrayfun(@(t) mean(metrics_mz.Aggregate.(metricName)(metrics_mz.Aggregate.trialInBlock == t), 'omitnan'), uniqueTrials);
        %     sems = arrayfun(@(t) calcSEM(metrics_mz.Aggregate.(metricName)(metrics_mz.Aggregate.trialInBlock == t)), uniqueTrials);
        %     errorbar(uniqueTrials, means, sems, 'o-', 'Color', mazeColors{3}, 'DisplayName', 'Aggregate', 'MarkerFaceColor', mazeColors{3});
        %
        %     title(metricName);
        %     xlabel('Trial in Block');
        %     ylabel(metricName);
        %     legend show;
        %
        %     hold off;
        % end
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


end

disp('done, return'), return



function metrics_mz = ExtractMetrics(mazeData, mazeType)

% mazeData will contain the data from Maze1 or Maze1Repeat for each session

if do_AcrossSessionAnalysis_SortedByMazeType
    % If sorting/filtering by Maze Type, extract metrics accordingly.
    % Assuming mazeTurnsLength first element can be used to distinguish types.
    if mazeData.mazeTurnsLength(j,1) == someConditionForMazeType
        metrics_mz = AppendMetrics(metrics_mz, mazeType, mazeData, j);
    end
    % Additional logic to handle different maze types can be added here.

elseif do_AcrossSessionAnalysis_SortedBySession
    % FILTERED BY RES INDEX
    % If sorting/filtering by Session, you might need session-specific data
    % which should be included in mazeData or handled externally.

elseif do_AcrossSessionAnalysis_SortedByWeekday
    % NEED TO FIND ALL SESSION WEEKDAYS {''}
    % Similar to session, handling based on weekday.

elseif (do_AcrossSessionAnalysis_LearningCurve || do_WithinSessionAnalysis)
    % FILTERED BY THE TRIAL COUNT IN BLOCK
    % If handling learning curves or within session analysis,
    % it could be direct extraction without specific conditions.
    metrics_mz = AppendMetrics(metrics_mz, mazeType, mazeData, j);
end
end


function metrics_mz = AppendMetrics(metrics_mz, mazeType, mazeData, metricsIndex, mazeDataIndex)

metrics_mz.(mazeType)   .normalizedTotalErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedTotalErrors{metricsIndex}, mazeData.totalErrors(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];
metrics_mz.(mazeType).normalizedDurations{metricsIndex} = [metrics_mz.(mazeType).normalizedDurations{metricsIndex}, mazeData.mazeDuration(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];
metrics_mz.(mazeType).normalizedRuleBreakingErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedRuleBreakingErrors{metricsIndex}, mazeData.ruleBreakingErrors(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];
metrics_mz.(mazeType).normalizedRuleAbidingErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedRuleAbidingErrors{metricsIndex}, mazeData.ruleAbidingErrors(mazeDataIndex) ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];

sumOfPerseverativeErrors = mazeData.perseverativeRuleBreakingErrors(mazeDataIndex) + ...
    mazeData.perseverativeRuleAbidingErrors(mazeDataIndex) + ...
    mazeData.perseverativeRetouchErroneous(mazeDataIndex);

metrics_mz.(mazeType).normalizedPerseverativeErrors{metricsIndex} = [metrics_mz.(mazeType).normalizedPerseverativeErrors{metricsIndex}, sumOfPerseverativeErrors ./ mazeData.mazeTurnsLength(mazeDataIndex,2)'];

end

function metrics_mz = AggregateMetrics(metrics_mz)
fields = fieldnames(metrics_mz.Maze1);
for f = 1:length(fields)
    fieldName = fields{f};
    metrics_mz.Aggregate.(fieldName) = [metrics_mz.Maze1.(fieldName), metrics_mz.Maze1Repeat.(fieldName)];
end
end


