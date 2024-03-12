
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


do_LoadData = 0;
do_AnalysisBasic = 1;
do_AnalysisBasicPlot = 1;
do_AcrossSessionAnalysis_SortedByMazeType = 0;
do_AcrossSessionAnalysis_SortedBySession = 0;
do_WithinSessionAnalysis = 1;

DO_SAVEFIGURES = 1;
PLOT_ACROSSSESSIONS_SortedByMazeType = 0;
PLOT_ACROSSSESSIONS_SortedBySession = 0;
PLOT_INDIVIDUALSESSIONS = 1;
% --- --- --- --- --- --- --- ---
% --- folder with preprocecced data
% --- --- --- --- --- --- --- ---
HOME_FOLDER = ['/Users/seema/Desktop/m-use/M-USE_Analysis' filesep];
%'/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS'
MUSEMATFOLDERNames = {}; iResultFolder = ''; iResultFile = '';


%SessionID_MZ = {'Frey_MZ_all_01'};
SessionID_MZ = {'Wotan_MZ_all_01'};

% MUSEMATFOLDERNames = 'MUSEMAT01_VS_EC_FL_MZG_Frey'
switch SessionID_MZ{1}
    case {'Frey_MZ_all_01'}
        MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_FL_MZG_Frey'
        % MUSEMATFOLDERNames{2} = 'MUSEMAT01_VS_EC_FL_MRT_Frey'
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Frey_01_20240208'; %'WM_01_20231205';%'WM_01_20230907';
    case {'Wotan_MZ_all_01'}
        MUSEMATFOLDERNames{1} =  'MUSEMAT01_WM_EC_FL_MZG_Wotan';
        %MUSEMATFOLDERNames{2} =  'MUSEMAT01_VS_EC_FL_MRT_Wotan';
        iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
        iResultFile   = 'MZ_Wotan_01_20240304';%'WM_01_20230907';
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
        % Identify unique mazes by turns and length
        mazeTurnLength = vertcat(res.mazeTurnsLength);
        uniqueMazes = unique(mazeTurnLength, 'rows');

        % Initialize metrics structure
        metrics_mz = struct();
        metrics_mz.subjectName = res(1).subjectName{1};
        metrics_mz.numTurnsTiles = uniqueMazes;
        metrics_mz.num = zeros(size(uniqueMazes, 1), 1);
        metrics_mz.mazeDefNames = cell(size(uniqueMazes, 1), 1);
        metrics_mz.sessionNumBlockNum = cell(size(uniqueMazes, 1), 1);
        metrics_mz.durations = cell(size(uniqueMazes, 1), 1);
        metrics_mz.completed = cell(size(uniqueMazes, 1), 1);
        metrics_mz.totalErrors = cell(size(uniqueMazes, 1), 1);
        metrics_mz.ruleBreakingErrors = cell(size(uniqueMazes, 1), 1);
        metrics_mz.ruleAbidingErrors = cell(size(uniqueMazes, 1), 1);
        metrics_mz.perseverationErrors = cell(size(uniqueMazes, 1), 1);
        metrics_mz.retouchErrors = cell(size(uniqueMazes, 1), 1);
        metrics_mz.retouchCorrect = cell(size(uniqueMazes, 1), 1);

        % Analyze each unique maze configuration
        for iL = 1:size(uniqueMazes, 1)
            sessionID = [];
            sessionIDTrials = {};
            iMazeDefNames = {};

            % Collect sessions and blocks for the current maze configuration
            for iDD = 1:length(res)
                if isempty(res(iDD).mazeTurnsLength), continue; end
                sel = find(res(iDD).mazeTurnsLength(:,1) == uniqueMazes(iL,1) & res(iDD).mazeTurnsLength(:,2) == uniqueMazes(iL,2));
                if ~any(sel), continue; end
                sessionID(end+1) = iDD;
                sessionIDTrials{end+1} = sel;
            end

            if isempty(sessionID), continue; end

            metrics_mz.num(iL) = numel(sessionID);

            % Temporary storage for metrics to be aggregated
            durations = [];
            completed = [];
            totalErrors = [];
            ruleAbidingErrors = [];
            ruleBreakingErrors = [];
            perseverationErrors = [];

            % Iterate through selected sessions and blocks to collect metrics
            for iSession = 1:numel(sessionID)
                iD = sessionID(iSession);
                trials = sessionIDTrials{iSession};
                for iTrial = trials'
                    % Assuming mazeDefName, mazeDuration, and sliderBarFilled are arrays with one element per trial
                    iMazeDefNames{end+1} = res(iD).mazeDefName{iTrial};
                    durations = [durations, res(iD).mazeDuration(iTrial)];
                    completed = [completed, res(iD).sliderBarFilled(iTrial)];
                    totalErrors = [totalErrors, res(iD).totalErrors(iTrial)];
                    ruleAbidingErrors = [ruleAbidingErrors, res(iD).ruleAbidingErrors(iTrial)];
                    ruleBreakingErrors = [ruleBreakingErrors, res(iD).ruleBreakingErrors(iTrial)];

                    sumOfPerseverativeErrors = res(iD).perseverativeRuleBreakingErrors(iTrial) + ...
                        res(iD).perseverativeRuleAbidingErrors(iTrial) + ...
                        res(iD).perseverativeRetouchErroneous(iTrial);
                    perseverationErrors = [perseverationErrors, sumOfPerseverativeErrors];

                end
            end

            % Store aggregated metrics for the current maze configuration
            metrics_mz.mazeDefNames{iL} = iMazeDefNames;
            metrics_mz.durations{iL} = durations;
            metrics_mz.completed{iL} = completed;
            metrics_mz.totalErrors{iL} = totalErrors;
            metrics_mz.ruleAbidingErrors{iL} = ruleAbidingErrors;
            metrics_mz.ruleBreakingErrors{iL} = ruleBreakingErrors;
            metrics_mz.perseverationErrors{iL} = perseverationErrors;
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
                metrics_mz1(i).normalizedDurations = mean((res(i).data.Maze1.mazeDuration ./ pathLengthMaze1), 'omitnan');
                metrics_mz1(i).normalizedTotalErrors = mean((res(i).data.Maze1.totalErrors ./ pathLengthMaze1), 'omitnan');
                metrics_mz1(i).normalizedRuleAbidingErrors = mean((res(i).data.Maze1.ruleAbidingErrors ./ pathLengthMaze1), 'omitnan');
                metrics_mz1(i).normalizedRuleBreakingErrors = mean((res(i).data.Maze1.ruleBreakingErrors ./ pathLengthMaze1), 'omitnan');

                sumOfPerseverativeErrors = res(i).data.Maze1.perseverativeRuleBreakingErrors + ...
                    res(i).data.Maze1.perseverativeRuleAbidingErrors + ...
                    res(i).data.Maze1.perseverativeRetouchErroneous;
                metrics_mz1(i).normalizedPerseverationErrors = mean((sumOfPerseverativeErrors ./ pathLengthMaze1), 'omitnan');
            end

            if isfield(res(i).data, 'Maze1Repeat')
                pathLengthMaze1Repeat = res(i).data.Maze1Repeat.mazeTurnsLength(:,2)';
                metrics_mz1Repeat(i).normalizedDurations = mean((res(i).data.Maze1Repeat.mazeDuration ./ pathLengthMaze1Repeat), 'omitnan');
                metrics_mz1Repeat(i).normalizedTotalErrors = mean((res(i).data.Maze1Repeat.totalErrors ./ pathLengthMaze1Repeat), 'omitnan');
                metrics_mz1Repeat(i).normalizedRuleAbidingErrors = mean((res(i).data.Maze1Repeat.ruleAbidingErrors ./ pathLengthMaze1Repeat), 'omitnan');
                metrics_mz1Repeat(i).normalizedRuleBreakingErrors = mean((res(i).data.Maze1Repeat.ruleBreakingErrors ./ pathLengthMaze1Repeat), 'omitnan');

                sumOfPerseverativeErrors = res(i).data.Maze1Repeat.perseverativeRuleBreakingErrors + ...
                    res(i).data.Maze1Repeat.perseverativeRuleAbidingErrors + ...
                    res(i).data.Maze1Repeat.perseverativeRetouchErroneous;
                metrics_mz1Repeat(i).normalizedPerseverationErrors = mean((sumOfPerseverativeErrors ./ pathLengthMaze1Repeat), 'omitnan');
            end
        end
    end

    if do_WithinSessionAnalysis == 1
        % Initialize metrics structure
        metrics_mz1 = res(6).data.Maze1;
        metrics_mz1Repeat = res(6).data.Maze1Repeat;
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


    if PLOT_ACROSSSESSIONS_SortedByMazeType==1

        iSubject = metrics_mz.subjectName;
        metrics_mz.numTurnsTiles
        minmaxTurns = [3 5];
        minmazTiles = [14 16]; iMazeLabel = 'T3.5/L14.16'
        %minmazTiles = [10 12]; iMazeLabel = 'T2/L10.12'
        sel =  find(metrics_mz.numTurnsTiles(:,1) >=minmaxTurns(1) & metrics_mz.numTurnsTiles(:,1) <=minmaxTurns(2)...
            & metrics_mz.numTurnsTiles(:,2) >=minmazTiles(1) & metrics_mz.numTurnsTiles(:,2) <=minmazTiles(2));

        % --- n mazes of this complextiy type:
        iNMazes = sum(metrics_mz.num(sel,1));
        % Find the maximum duration length for padding
        maxDurationLength = max(cellfun(@length, metrics_mz.durations(sel)));

        % Initialize structure for normalized metrics
        normalizedMetrics = struct();
        normalizedMetrics.durations = cell(1, length(sel));
        normalizedMetrics.totalErrors = cell(1, length(sel));
        normalizedMetrics.ruleBreakingErrors = cell(1, length(sel));
        normalizedMetrics.ruleAbidingErrors = cell(1, length(sel));
        normalizedMetrics.perseverationErrors = cell(1, length(sel));


        for j = 1:length(sel)
            nTiles = metrics_mz.numTurnsTiles(sel(j), 2);
            durations = metrics_mz.durations{sel(j)};
            totalErrors = metrics_mz.totalErrors{sel(j)};
            ruleBreakingErrors = metrics_mz.ruleBreakingErrors{sel(j)};
            ruleAbidingErrors = metrics_mz.ruleAbidingErrors{sel(j)};
            perseverationErrors = metrics_mz.perseverationErrors{sel(j)};

            % Normalize and store the metrics
            normalizedMetrics.durations{j} = durations ./ nTiles;
            normalizedMetrics.totalErrors{j} = totalErrors ./ nTiles;
            normalizedMetrics.ruleBreakingErrors{j} = ruleBreakingErrors ./ nTiles;
            normalizedMetrics.ruleAbidingErrors{j} = ruleAbidingErrors ./ nTiles;
            normalizedMetrics.perseverationErrors{j} = perseverationErrors ./ nTiles;
        end

        figure('Color', 'w');
        set(gcf, 'Position', [100, 100, 1800, 600]); % Adjust figure size to fit 3 plots in a row
        nc = 3; % Number of columns in subplot grid
        nr = 1; % Number of rows in subplot grid

        % Preparing mazeConditions labels for x-axis
        mazeConditions = strings(1, length(sel));
        for j = 1:length(sel)
            mazeConditions(j) = "T" + num2str(metrics_mz.numTurnsTiles(sel(j),1)) + ".L" + num2str(metrics_mz.numTurnsTiles(sel(j),2));
        end

        % Function to calculate SEM
        calcSEM = @(data) std(data) / sqrt(length(data));

        % Plot 1: Durations Plot
        subplot(nr, nc, 1);
        hold on;
        for j = 1:length(sel)
            durations = cell2mat(normalizedMetrics.durations(j));
            medianDur = median(durations);
            semDur = calcSEM(durations);  % Calculate SEM
            errorbar(j, medianDur, semDur, 'o');
        end
        xlim([0.5, length(sel) + 0.5]);
        ylabel('Normalized Duration');
        xlabel('Condition');
        title('Normalized Durations Across Conditions');
        set(gca, 'XTick', 1:length(sel), 'XTickLabel', mazeConditions);
        hold off;

        % Plot 2: Rule-Breaking and Rule-Abiding Errors overlaid
        subplot(nr, nc, 2);
        hold on;
        for j = 1:length(sel)
            % Rule-Breaking Errors
            ruleBreakingErrors = cell2mat(normalizedMetrics.ruleBreakingErrors(j));
            medianRB = median(ruleBreakingErrors);
            semRB = calcSEM(ruleBreakingErrors);  % Calculate SEM

            % Rule-Abiding Errors
            ruleAbidingErrors = cell2mat(normalizedMetrics.ruleAbidingErrors(j));
            medianRA = median(ruleAbidingErrors);
            semRA = calcSEM(ruleAbidingErrors);  % Calculate SEM

            % Plotting with SEM
            errorbar(j-0.1, medianRB, semRB, 'o', 'MarkerFaceColor', 'r', 'Color', 'r'); % Slightly offset for visibility
            errorbar(j+0.1, medianRA, semRA, 'o', 'MarkerFaceColor', 'b', 'Color', 'b'); % Slightly offset for visibility
        end
        xlim([0.5, length(sel) + 0.5]);
        ylabel('Errors');
        xlabel('Condition');
        title('Rule-Breaking vs. Rule-Abiding Errors');
        legend('Rule-Breaking', 'Rule-Abiding', 'Location', 'best');
        set(gca, 'XTick', 1:length(sel), 'XTickLabel', mazeConditions);
        hold off;

        % Plot 3: Total Errors
        subplot(nr, nc, 3);
        hold on;
        for j = 1:length(sel)
            totalErrors = cell2mat(normalizedMetrics.totalErrors(j));
            medianTE = median(totalErrors);
            semTE = calcSEM(totalErrors);  % Calculate SEM
            errorbar(j, medianTE, semTE, 'o', 'MarkerFaceColor', 'k', 'Color', 'k');
        end
        xlim([0.5, length(sel) + 0.5]);
        ylabel('Normalized Total Errors');
        xlabel('Condition');
        title('Total Errors Across Conditions');
        set(gca, 'XTick', 1:length(sel), 'XTickLabel', mazeConditions);
        hold off;
    end

   if PLOT_ACROSSSESSIONS_SortedBySession == 1
    % Extract session dates from dataset names
    sessionDates = cellfun(@(x) regexp(x, 'Session_(\d+_\d+)', 'tokens'), {res.dataset}, 'UniformOutput', false);
    sessionDates = cellfun(@(x) x{1}{1}, sessionDates, 'UniformOutput', false); % Convert tokens to strings
    sessionDatesFormatted = datetime(sessionDates, 'InputFormat', 'MM_dd', 'Format', 'MM/dd'); % Convert to datetime

    % Initialize arrays to hold the means for plotting
    normalizedDurationsMean = [];
    normalizedTotalErrorsMean = [];
    normalizedRuleAbidingErrorsMean = [];
    normalizedRuleBreakingErrorsMean = [];

    % Collect the means from the metrics_mz1 and metrics_mz1Repeat structures
    for i = 1:length(metrics_mz1)
        normalizedDurationsMean = [normalizedDurationsMean; metrics_mz1(i).normalizedDurations, metrics_mz1Repeat(i).normalizedDurations];
        normalizedTotalErrorsMean = [normalizedTotalErrorsMean; metrics_mz1(i).normalizedTotalErrors, metrics_mz1Repeat(i).normalizedTotalErrors];
        normalizedRuleAbidingErrorsMean = [normalizedRuleAbidingErrorsMean; metrics_mz1(i).normalizedRuleAbidingErrors, metrics_mz1Repeat(i).normalizedRuleAbidingErrors];
        normalizedRuleBreakingErrorsMean = [normalizedRuleBreakingErrorsMean; metrics_mz1(i).normalizedRuleBreakingErrors, metrics_mz1Repeat(i).normalizedRuleBreakingErrors];
    end

    % Create a figure for the plots
    figure;

    % Specify dates for vertical lines
    specificDates = ["02/12", "02/19", "02/26", "03/04"];
    year = 2024; % Assuming the year is 2024, adjust accordingly
    datesForLines = datetime(specificDates + "/" + num2str(year), 'InputFormat', 'MM/dd/yyyy');

    % Loop through each subplot to plot data
    for plotIndex = 1:4
        subplot(2, 2, plotIndex);
        hold on;

        % Select the appropriate dataset for each subplot
        switch plotIndex
            case 1
                dataToPlot = normalizedDurationsMean;
                plotTitle = 'Normalized Durations Mean';
            case 2
                dataToPlot = normalizedTotalErrorsMean;
                plotTitle = 'Normalized Total Errors Mean';
            case 3
                dataToPlot = normalizedRuleAbidingErrorsMean;
                plotTitle = 'Normalized Rule Abiding Errors Mean';
            case 4
                dataToPlot = normalizedRuleBreakingErrorsMean;
                plotTitle = 'Normalized Rule Breaking Errors Mean';
        end

        % Plot lines and scatter for the selected dataset
        for seriesIndex = 1:size(dataToPlot, 2)
            if seriesIndex == 1
                displayName = 'Maze1';
            else
                displayName = 'Maze1Repeat';
            end
            plot(sessionDatesFormatted, dataToPlot(:, seriesIndex), 'LineWidth', 1.5, 'DisplayName', displayName);
            scatter(sessionDatesFormatted, dataToPlot(:, seriesIndex), 'filled', 'HandleVisibility', 'off');
        end

        % Draw vertical lines at specific dates
        for xValue = datesForLines
            xline(xValue, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
        end

        title(plotTitle);
        xlabel('Session Date');
        ylabel('Value');
        xticks(sessionDatesFormatted);  % Set x-axis ticks at each date
        xtickformat('MM/dd');  % Format the tick labels
        xtickangle(45);  % Rotate labels for readability
        xlim([min(sessionDatesFormatted) max(sessionDatesFormatted)]);  % Adjust x-axis limits

        % Add legend in the first plot only
        if plotIndex == 1
            legend('Location', 'best');
        end

        hold off;
    end

    % Adjust figure size and display
    set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust size as needed
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
