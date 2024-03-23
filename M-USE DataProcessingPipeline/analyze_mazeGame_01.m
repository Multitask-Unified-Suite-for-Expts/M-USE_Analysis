
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
do_WithinSessionAnalysis = 1;

DO_SAVEFIGURES = 0;
PLOT_ACROSSSESSIONS_SortedByMazeType = 0;
PLOT_ACROSSSESSIONS_SortedBySession = 0;
PLOT_INDIVIDUALSESSIONS = 1;
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
                    if ~isfield(res(iD).data, strrep(label{1}, '_', ''))
                        res(iD).data.(strrep(label{1}, '_', '')) = struct();
                    end

                    trialData = in.dat.trialData{idx};
                    blockData = in.dat.blockData{idx};
                    blockDef = table2struct(in.dat.cfg_BlockDef{idx});
                    iT = 0;
                    for jj = 1:length(trialData)
                        if trialData(jj).AbortCode ~= 0
                            continue; % Skip to the next iteration of the loop
                        end
                        % Assign trial data, block data, and block definitions
                        currentData = res(iD).data.(strrep(label{1}, '_', ''));

                        iT = iT + 1;
                        currentData.subjectName{iT} = trialData(jj).SubjectID;%  trialData(iT).SubjectNum;

                        mazeStruct = jsondecode(blockDef(trialData(jj).BlockCount).MazeDef);

                        currentData.mazeDefName{iT} = mazeStruct.mName;
                        nTurns = mazeStruct.mNumTurns;
                        mLength = mazeStruct.mNumSquares;
                        currentData.mazeTurnsLength(iT,1:2) = [nTurns, mLength];
                        currentData.blockNum(iT) = trialData(jj).BlockCount;
                        currentData.trialInBlock(iT) = trialData(jj).TrialCount_InBlock;

                        currentData.mazeDuration(iT) = trialData(jj).MazeDuration;
                        currentData.sliderBarFilled(iT) = strcmp(trialData(jj).SliderBarFilled,'True');
                        currentData.totalErrors(iT) = trialData(jj).TotalErrors;
                        currentData.selectedTiles{iT} = strsplit(trialData(jj).SelectedTiles, ',');
                        currentData.correctTouches(iT) = trialData(jj).CorrectTouches;
                        currentData.retouchCorrect(iT) = trialData(jj).RetouchCorrect;
                        currentData.retouchErroneous(iT) = trialData(jj).RetouchErroneous;
                        currentData.backTrackingErrors(iT) = trialData(jj).BacktrackingErrors;
                        currentData.ruleAbidingErrors(iT) = trialData(jj).Rule_AbidingErrors;
                        currentData.ruleBreakingErrors(iT) = trialData(jj).Rule_BreakingErrors;
                        currentData.perseverativeRetouchErroneous(iT) = trialData(jj).PerseverativeRetouchErrors;
                        currentData.perseverativeBackTrackingErrors(iT) = trialData(jj).PerseverativeBackTrackErrors;
                        currentData.perseverativeRuleAbidingErrors(iT) = trialData(jj).PerseverativeRuleAbidingErrors;
                        currentData.perseverativeRuleBreakingErrors(iT) = trialData(jj).PerseverativeRuleBreakingErrors;
                        % Save the collected data back to the main structure
                        res(iD).data.(strrep(label{1}, '_', '')) = currentData;
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
% --- --- --- --- --- --- --- --- --- --- --- --- ---
if do_AnalysisBasic == 1

    % Load results if they are not yet loaded
    if ~exist('res', 'var')
        load(fullfile(RESULTFOLDER, iResultFile))
    end

    if do_AcrossSessionAnalysis_SortedByMazeType == 1
        % Extracting unique mazes based on new structure
        allMazes = [];
        for i = 1:length(res)
            if isfield(res(i).data, 'Maze1')
                allMazes = [allMazes; res(i).data.Maze1.mazeTurnsLength];
            end
            if isfield(res(i).data, 'Maze1Repeat')
                allMazes = [allMazes; res(i).data.Maze1Repeat.mazeTurnsLength];
            end
        end

        uniqueMazes = unique(allMazes, 'rows');

        % Filter uniqueMazes to include only those with a 16 length
        % uniqueMazes= uniqueMazes(uniqueMazes(:,2) == 16, :);

        % Initialize metrics structure
        metrics_mz = struct();
        metrics_mz.subjectName = res(1).data.Maze1.subjectName{1};
        metrics_mz.numTurnsTiles = uniqueMazes;
        metrics_mz.num = zeros(size(uniqueMazes, 1), 1);
        metrics_mz.mazeDefNames = cell(size(uniqueMazes, 1), 1);
        metrics_mz.sessionNumBlockNum = cell(size(uniqueMazes, 1), 1);
        % New fields for normalized metrics
        metrics_mz.normalizedDurations = cell(size(uniqueMazes, 1), 3); % First, Repeat, Aggregate
        metrics_mz.normalizedTotalErrors = cell(size(uniqueMazes, 1), 3); % and so on
        metrics_mz.normalizedRuleAbidingErrors = cell(size(uniqueMazes, 1), 3);
        metrics_mz.normalizedRuleBreakingErrors = cell(size(uniqueMazes, 1), 3);
        metrics_mz.normalizedPerseverationErrors = cell(size(uniqueMazes, 1), 3);

        % Analyze each unique maze configuration
        for iL = 1:size(uniqueMazes, 1)
            % Variables to accumulate metrics
            durationsFirst = [];
            durationsRepeat = [];
            totalErrorsFirst = [];
            totalErrorsRepeat = [];
            % Similar variables for other metrics...

            for iD = 1:length(res)
                % Handle Maze1
                if isfield(res(iD).data, 'Maze1')
                    maze1 = res(iD).data.Maze1;
                    sel = find(maze1.mazeTurnsLength(:,1) ==  uniqueMazes(iL,1) & maze1.mazeTurnsLength(:,2) == uniqueMazes(iL,2));
                    if ~isempty(sel)
                        durationsFirst = [durationsFirst, maze1.mazeDuration(sel) ./ maze1.mazeTurnsLength(sel,2)'];
                        totalErrorsFirst = [totalErrorsFirst, maze1.totalErrors(sel) ./ maze1.mazeTurnsLength(sel,2)'];                    end
                end

                % Handle Maze1Repeat
                if isfield(res(iD).data, 'Maze1Repeat')
                    maze1Repeat = res(iD).data.Maze1Repeat;
                    sel = find(maze1Repeat.mazeTurnsLength(:,1) ==  uniqueMazes(iL,1) & maze1Repeat.mazeTurnsLength(:,2) == uniqueMazes(iL,2));
                    if ~isempty(sel)
                        durationsRepeat = [durationsRepeat, maze1Repeat.mazeDuration(sel) ./ maze1Repeat.mazeTurnsLength(sel,2)'];
                        totalErrorsRepeat = [totalErrorsRepeat, maze1Repeat.totalErrors(sel) ./ maze1Repeat.mazeTurnsLength(sel,2)'];
                    end
                end
            end

            % Aggregate metrics (First, Repeat, and Combined)
            metrics_mz.normalizedDurations{iL, 1} = durationsFirst;
            metrics_mz.normalizedDurations{iL, 2} = durationsRepeat;
            metrics_mz.normalizedDurations{iL, 3} = [durationsFirst, durationsRepeat];

            metrics_mz.normalizedTotalErrors{iL, 1} = totalErrorsFirst;
            metrics_mz.normalizedTotalErrors{iL, 2} = totalErrorsRepeat;
            metrics_mz.normalizedTotalErrors{iL, 3} = [totalErrorsFirst, totalErrorsRepeat];

        end
        % Save analysis results
        save(fullfile(RESULTFOLDER, iResultFileMetrics), 'res', 'metrics_mz', '-v7.3')
        disp(['Saved metrics mat-file of n=', num2str(length(res)), ' datasets in ', iResultFileMetrics]);


    end
    if do_AcrossSessionAnalysis_SortedBySession == 1
        % Initialize structures to store the aggregated metrics
        metrics_mz1 = struct('normalizedDurations', {}, 'normalizedTotalErrors', {}, 'normalizedRuleAbidingErrors', {}, 'normalizedRuleBreakingErrors', {}, 'normalizedPerseverationErrors', {});
        metrics_mz1Repeat = struct('normalizedDurations', {}, 'normalizedTotalErrors', {}, 'normalizedRuleAbidingErrors', {}, 'normalizedRuleBreakingErrors', {}, 'normalizedPerseverationErrors', {});

        for i = 1:length(res)
            if isfield(res(i).data, 'Maze1')
                pathLengthMaze1 = res(i).data.Maze1.mazeTurnsLength(:,2)';
                metrics_mz1(i).meanNormalizedDurations = mean((res(i).data.Maze1.mazeDuration ./ pathLengthMaze1), 'omitnan');
                metrics_mz1(i).meanNormalizedTotalErrors = mean((res(i).data.Maze1.totalErrors ./ pathLengthMaze1), 'omitnan');
                metrics_mz1(i).meanNormalizedRuleAbidingErrors = mean((res(i).data.Maze1.ruleAbidingErrors ./ pathLengthMaze1), 'omitnan');
                metrics_mz1(i).meanNormalizedRuleBreakingErrors = mean((res(i).data.Maze1.ruleBreakingErrors ./ pathLengthMaze1), 'omitnan');

                sumOfPerseverativeErrors = res(i).data.Maze1.perseverativeRuleBreakingErrors + ...
                    res(i).data.Maze1.perseverativeRuleAbidingErrors + ...
                    res(i).data.Maze1.perseverativeRetouchErroneous;
                metrics_mz1(i).meanNormalizedPerseverationErrors = mean((sumOfPerseverativeErrors ./ pathLengthMaze1), 'omitnan');
            end

            if isfield(res(i).data, 'Maze1Repeat')
                pathLengthMaze1Repeat = res(i).data.Maze1Repeat.mazeTurnsLength(:,2)';
                metrics_mz1Repeat(i).meanNormalizedDurations = mean((res(i).data.Maze1Repeat.mazeDuration ./ pathLengthMaze1Repeat), 'omitnan');
                metrics_mz1Repeat(i).meanNormalizedTotalErrors = mean((res(i).data.Maze1Repeat.totalErrors ./ pathLengthMaze1Repeat), 'omitnan');
                metrics_mz1Repeat(i).meanNormalizedRuleAbidingErrors = mean((res(i).data.Maze1Repeat.ruleAbidingErrors ./ pathLengthMaze1Repeat), 'omitnan');
                metrics_mz1Repeat(i).meanNormalizedRuleBreakingErrors = mean((res(i).data.Maze1Repeat.ruleBreakingErrors ./ pathLengthMaze1Repeat), 'omitnan');

                sumOfPerseverativeErrors = res(i).data.Maze1Repeat.perseverativeRuleBreakingErrors + ...
                    res(i).data.Maze1Repeat.perseverativeRuleAbidingErrors + ...
                    res(i).data.Maze1Repeat.perseverativeRetouchErroneous;
                metrics_mz1Repeat(i).meanNormalizedPerseverationErrors = mean((sumOfPerseverativeErrors ./ pathLengthMaze1Repeat), 'omitnan');
            end
        end
    end

    if do_WithinSessionAnalysis == 1
        % Initialize metrics structure
        metrics_mz1 = res(13).data.Maze1;
        metrics_mz1Repeat = res(13).data.Maze1Repeat;
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


    if PLOT_ACROSSSESSIONS_SortedByMazeType == 1

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


    if PLOT_ACROSSSESSIONS_SortedBySession == 1
        % Extract session dates from dataset names
        sessionDates = cellfun(@(x) regexp(x, 'Session_(\d+_\d+)', 'tokens'), {res.dataset}, 'UniformOutput', false);
        sessionDates = cellfun(@(x) x{1}{1}, sessionDates, 'UniformOutput', false); % Convert tokens to strings
        sessionDatesFormatted = datetime(sessionDates, 'InputFormat', 'MM_dd', 'Format', 'MM/dd'); % Convert to datetime

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


    if(PLOT_INDIVIDUALSESSIONS)
        % Prepare a single figure with two subplots
        figure;
        % Extract relevant fields from both structs
        blockNumsInitial = metrics_mz1.blockNum;
        trialInBlockInitial = metrics_mz1.trialInBlock;
        totalErrorsInitial = metrics_mz1.totalErrors;
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
        ylabel('Normalized Maze Duration');
        title('Normalized Maze Durations');
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
