

%
% This script contains some basic analysis of the Working Memory (WM) task.
% It assumes that the script analyze_preprocessSessions_01.m was run to
% save trialdata from 1 or more sessions containing the VS task

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
%    M-USE matlab scripts are distributed in the hope that they will be
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
%    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    A copy of the GNU General Public License can be found
%    at http://www.gnu.org/licenses/
%

% id 2023-08-30 analyze_workingMemory_01
% FOLDER_DATA = '/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/#ms_CognitiveGames_Kiosk/GAMES_MATLAB_ANALYSIS/MUSE_analysis_Maze/DATA/'

do_LoadData = 1;
do_AnalysisBasic = 1
do_AnalysisBasicPlot = 1;

% --- --- --- --- --- --- --- ---
% --- folder with preprocecced data
% --- --- --- --- --- --- --- ---
HOME_FOLDER = ['/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS' filesep];
                %'/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS'
MUSEMATFOLDERNames = {}; iResultFolder = ''; iResultFile = '';


% SessionID_WM = {'Wotan_WM_all_01'}
SessionID_WM = {'Frey_WM_all_01'}
%SessionID_WM = {'Wotan_WM_all_01'}


% MUSEMATFOLDERNames = 'MUSEMAT01_VS_EC_FL_MZG_Frey'
if strcmp(SessionID_WM{1},'Frey_WM_all_01')
    MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_FL_MZG_Frey'
    
    iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
    iResultFile   = 'WM_Frey_01_20240129';%'WM_01_20231205';%'WM_01_20230907';
end
if strcmp(SessionID_WM{1},'Wotan_WM_all_01')
    MUSEMATFOLDERNames{1} =  'MUSEMAT01_WM_EC_FL_MZG_Wotan';
    % MUSEMATFOLDERNames{2} =  'MUSEMAT01_VS_EC_FL_MRT_Wotan';
    iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
    iResultFile   = 'WM_Wotan_01_20231205';%'WM_01_20230907';
end

iMonkey = SessionID_WM{1}(1:findstr(SessionID_WM{1},'_')-1);

% --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
% --- The files generates by this script and the folder name for them:
iResultFileMetrics = [iResultFile '_Metrics'];
RESULTFOLDER = [HOME_FOLDER filesep iResultFolder '_MAT']; if ~exist(RESULTFOLDER), mkdir(RESULTFOLDER),end
METRICSFOLDER  = [RESULTFOLDER '_METRICS']; if ~exist(METRICSFOLDER), mkdir(METRICSFOLDER),end
FIGURE_Folder = [RESULTFOLDER '_FIG'];    if ~exist(FIGURE_Folder), mkdir(FIGURE_Folder), end





% --- --- ---
% --- Loading data and organizing them in a mat file
% --- --- ---
if do_LoadData==1

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

    res = [];
    res.dataset = [];
    res.datafolder = [];
    res.SubjectID =  [];
    res.nTrials  = [];

    % --- --- --- --- --- --- --- --- --- --- --- ---
    % --- specify which features to expect:
    % --- --- --- --- --- --- --- --- --- --- --- ---
    res.dimensionNums = {};
    res.dimensionNames = {};
    res.dimensionNums{1}=1:9;
    res.dimensionNames{1}={'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09'};
    res.dimensionNums{2}=10:18;
    res.dimensionNames{2}={'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09'};
    res.dimensionNums{3}=19:26;
    res.dimensionNames{3}{1}={'C6070014_6070014', 'C6070059_6070059', 'C6070106_6070106', 'C6070148_6070148', 'C6070194_6070194', 'C6070239_6070239', 'C6070287_6070287', 'C6070335_6070335'};
    res.dimensionNames{3}{2}={'C7070014_5000000', 'C7070059_5000000', 'C7070106_5000000', 'C7070148_5000000', 'C7070194_5000000', 'C7070240_5000000', 'C7070286_5000000', 'C7070335_5000000'};
    res.dimensionNums{4}=27:37;
    res.dimensionNames{4}={'A00_E01', 'A00_E02', 'A00_E03', 'A01_E00', 'A01_E01', 'A01_E02', 'A01_E03', 'A02_E00', 'A02_E01', 'A02_E02', 'A02_E03'};
    res.ObjectFeatureLabels = [res.dimensionNames{1} res.dimensionNames{2} res.dimensionNames{3}{1} res.dimensionNames{4}];
    res.DimensionVector = [zeros(1,length(res.dimensionNums{1}))+1, zeros(1,length(res.dimensionNums{2}))+2, zeros(1,length(res.dimensionNums{3}))+3, zeros(1,length(res.dimensionNums{4}))+4];

    res.TrialInBlock = [];
    res.TrialConditionLabel = [];
    res.ContextID = [];
    res.ITI = [];
    res.trialID = {};
    res.delayDuration = [];
    res.displaySampleDuration = [];

    res.BlockNum = [];
    res.TrialInExperiment = [];
    res.TrialInBlock = [];
    res.Outcome  = [];
    res.AllTokensCompleted  = [];
    res.TokenBarValue  = [];
    res.Accuracy = [];
    res.ReactionTime = [];
    res.AbortCode = [];

    res.chosenObjectVector = [];
    res.targetObjectVector = [];
    res.distractorVectors = [];

    res.chosenObjectPositionVector = [];
    res.targetObjectPositionVector = [];
    res.distractorPositionVectors = [];

    res.nDistractors = [];
    res.TDsharedFeatures = [];

    res.cnd_DelaySec = [];
    res.cnd_DistractorSetSize = [];

    for iD=1:length(datasets)

               in = load([datasetsDIR{iD} filesep datasets{iD}]); disp('.')
            
               % disp('hhh'), return
               if isempty(in.dat.taskLabel), disp(sprintf('empty tasklabel for iD=%d, %s',iD,datasets{iD})), continue, end

                % --- --- --- --- --- --- --- --- --- --- --- 
                % --- get all instances of the current task
                % --- --- --- --- --- --- --- --- --- --- --- 
        % idx_all = find(contains(in.dat.taskLabel,'FlexLearning'));
        % idx_all = find(contains(in.dat.taskLabel,'VisualSearch'));
        % idx_all = find(contains(in.dat.taskLabel,'EffortControl'));
        % idx_all = find(contains(in.dat.taskLabel,'MazeGuided'));
        % idx_all = find(contains(in.dat.taskLabel,'Maze'));
        idx_all = find(contains(in.dat.taskLabel,'WorkingMemory'));
        if isempty(idx_all), res(iD).log=[]; res(iD).log{end+1} = 'task isempty...'; continue, end
        for iIDX = 1:length(idx_all)
            idx = idx_all(iIDX);

% disp('hhh'), return
% ls /Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS/RES01_WM_EC_FL_MZG_VS_Wotan_MAT/
 %disp('hhh'), return
       disp(sprintf('processing iD=%d, %s',iD,datasets{iD})),
        if isempty(in.dat.trialData) & isempty(in.dat.blockData) ...
                | (length(in.dat.trialData) < idx)  | (length(in.dat.blockData) < idx)
            disp('EMPTY data'),
            continue, 
        end

        res(iD).datafolder = [datasetsDIR{iD}];
        res(iD).dataset = [datasets{iD}];
        res(iD).log = [];

        trialData = [];
        blockData = [];
        blockDef = [];
        stimDef = [];
        trialDef = [];
        if ~isempty(in.dat.trialData{idx})
            trialData = in.dat.trialData{idx};
        end

        if ~isempty(in.dat.blockData{idx})
            blockData = in.dat.blockData{idx};
        end

        if ~isempty(in.dat.cfg_BlockDef{idx})
            blockDef = table2struct(in.dat.cfg_BlockDef{idx});
        end
        if ~isempty(in.dat.cfg_StimDef{idx})
            stimDef  = table2struct(in.dat.cfg_StimDef{idx});
        end
        if ~isempty(in.dat.cfg_TrialDef{idx})
            trialDef = table2struct(in.dat.cfg_TrialDef{idx});
        end
        % dat.blockData{idx};
        % dat.frameData{idx};

        res(iD).nTrials  = length(trialData);

        % function FLData = get_trialData_FLUToken_02(iFileAndPath, iSessionName, trialdefname, trialData, blockData, blockjson, experimentID)
        for iT=1:res(iD).nTrials
            res(iD).SubjectID{iT}   = trialData(iT).SubjectID;
            res(iD).ITI(iT)         = trialData(iT).ITI_Duration;

            res(iD).trialID{iT}         = trialData(iT).TrialID;

            res(iD).delayDuration(iT)         = trialData(iT).Delay_Duration;
            res(iD).displaySampleDuration(iT)         = trialData(iT).DisplaySample_Duration;

            % TrialID: 'WMn1.TDSIM1_delay1.0' but Delay_Duration: 0.5171 and the
            % framedta suggests: (3542-3513)*(1000/60)
            % what is DisplaySample_Duration

            res(iD).BlockNum(iT)    = trialData(iT).BlockCount;
            if isnumeric(trialData(iT).TrialCount_InTask)
                res(iD).TrialInExperiment(iT) = trialData(iT).TrialCount_InTask;
                res(iD).TrialInBlock(iT)    = trialData(iT).TrialCount_InBlock;
            else
                error('TrialCount_InTask not numeric - why?');
            end

            res(iD).ContextID{iT}  = trialData(iT).ContextName; %Context
            res(iD).Outcome(iT) = trialData(iT).TokenChange;%tTokenRewardSize(1:nTrials);
            res(iD).AllTokensCompleted = strcmp(trialData(iT).TokenBarFull,'True');%TokenBarCompletedThisTrial
            res(iD).TokenBarValue = trialData(iT).TokenBarValue;

            res(iD).correctObjectPositionVector{iT} = [trialData(iT).SelectedLocation];
            res(iD).Accuracy(iT) = trialData(iT).CorrectSelection;
            res(iD).ReactionTime(iT) = trialData(iT).SearchDuration;
            res(iD).AbortCode(iT) = trialData(iT).AbortCode;

            % --- if the trial was aborted then continue
            if isnan(trialData(iT).SelectedStimCode),
                res(iD).chosenObjectVector(iT,:)            = nan(1,max(res(1).dimensionNums{end}));
                res(iD).chosenObjectPositionVector(iT,:)    = nan(1,3);
                res(iD).targetObjectVector(iT,:)            = nan(1,max(res(1).dimensionNums{end}));
                res(iD).targetObjectPositionVector(iT,:)    = nan(1,3);
                res(iD).distractorVectors{iT}               = [];
                res(iD).distractorPositionVectors{iT}       = [];
                continue,
            end


            % --- --- --- --- --- --- --- --- --- --- --- ---
            % --- for Quaddle 1.0 stimuli specified in the stimDef:
            % --- specify for each stimulus in this trial the features it carried
            % ---
            % --- requires: (1) stimulus indices in trialData : iStimIndices
            % ---           (2) the stimDef with index - name pairs : stimDef
            % ---           (3) dimension information about feature the filename
            % ---               refers to : dimensionNames, dimensionNums
            % --- --- --- --- --- --- --- --- --- --- --- ---
            iTrial = trialData(iT).TrialCount_InTask;

            % --- --- --- --- --- --- ---
            % --- specifc the chosen stimulus and chosen location
            % --- --- --- --- --- --- ---
            %iChosenStimIdx = trialData(iT).SelecteStimIndex; %SelectedStimIndex
            %NOTE: trialData(iT).SelectedStimCode is always 0 in working memory task....
            %        iChosenStimIdx = trialData(iT).SelectedStimCode; %SelectedStimIndex
            %        objectVectorInformation = get_objectVector_01(iChosenStimIdx, stimDef, res.dimensionNames, res.dimensionNums);
            %        res(iD).chosenObjectVector(iT,:) = objectVectorInformation{1};
            trialData(iT).SelectedLocation([findstr(trialData(iT).SelectedLocation,'(') findstr(trialData(iT).SelectedLocation,')')])=[];
            %        res(iD).chosenObjectPositionVector(iT,:) = str2num(trialData(iT).SelectedLocation);

            % --- --- --- --- --- --- ---
            % --- the first stimindex is the target (per definition)
            % --- --- --- --- --- --- ---

            %iStimIndices = str2num(trialDef(iTrial).TrialStimIndices);
            %iLocationsXYZ = str2num(trialDef(iTrial).TrialStimLocations);
            iStimIndices = str2num(trialDef(iTrial).SearchStimIndices);
            iLocationsXYZ = str2num(trialDef(iTrial).SearchStimLocations);

            % --- check chosen stimnulus was in stimulus list of stimDef file
            %if isempty(find(iChosenStimIdx==iStimIndices)), error('fatal error: not target indices (in trialData) found in trialDef stimindices'), end

            % --- the first stimindex is the target (per definition)
            iTargetStimIndex = iStimIndices(1);
            objectVectorInformation = get_objectVector_01(iTargetStimIndex, stimDef, res(1).dimensionNames, res(1).dimensionNums);
            res(iD).targetObjectVector(iT,:) = objectVectorInformation{1};
            res(iD).targetObjectPositionVector(iT,:) = iLocationsXYZ(1:3);

            % --- the second to last stimindices are distractors (per definition)
            res(iD).distractorVectors{iT} = [];
            res(iD).distractorPositionVectors{iT} = [];
            if length(iStimIndices)>1 % only find object vectors of distractors if there are distractors
                iStimIndices(1)=[];
                objectVectorInformation = get_objectVector_01(iStimIndices, stimDef, res(1).dimensionNames, res(1).dimensionNums);
                res(iD).distractorVectors{iT} = objectVectorInformation;
                iP=0;
                for p=1:3:length(iStimIndices)*3
                    iP=iP+1;
                    res(iD).distractorPositionVectors{iT}(iP,:) = iLocationsXYZ(p+[0 1 2]);
                end
                % --- if chosen object is this distractor its at the same location ...
                tmp = str2num(trialData(iT).SelectedLocation);
                if sum(tmp == res(iD).distractorPositionVectors{iT})==3
                    res(iD).chosenObjectVector(iT,:) = res(iD).distractorVectors{iT}{:};
                end

            end

            % --- if chosen object is target its at the same location ...
            tmp = str2num(trialData(iT).SelectedLocation);
            if sum(tmp == res(iD).targetObjectPositionVector)==3
                res(iD).chosenObjectVector(iT,:) = res(iD).targetObjectVector(iT,:);
            end



            % --- collect information about number of distractors and how many feature
            % --- were are shared between distractor and target ?

            res(iD).nDistractors(iT) = length(res(iD).distractorVectors{iT});
            res(iD).TDsharedFeatures(iT,1:4) = nan(1,4);
            if ~isempty(res(iD).distractorVectors{iT})
                tmp = [];
                for iDV=1:length(res(iD).distractorVectors{iT})
                    tmp(iDV) = sum(sum(cat(1,res(iD).targetObjectVector(iT,:), res(iD).distractorVectors{iT}{iDV}),1)==2);
                end
                res(iD).TDsharedFeatures(iT,1:length(tmp)) = tmp;

            end
            sprintf('done with dataset %d\n', iD);
        end

        end
        disp('done processing data\n');
    end
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    % --- save analysis results in RES file
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    save([RESULTFOLDER filesep iResultFile],'res','-V7.3')
    fprintf('saved analysis of n=%d datasets in %s \n',iD, iResultFile),
end



% dimensionNums: {[1 2 3 4 5 6 7 8 9]  [10 11 12 13 14 15 16 17 18]  [19 20 21 22 23 24 25 26]  [27 … ]}
% dimensionNames: {{1×9 cell}  {1×9 cell}  {1×8 cell}  {1×11 cell}}

% --- --- --- --- --- --- --- --- --- --- --- --- ---
% --- do some basic analysis
% --- --- --- --- --- --- --- --- --- --- --- --- ---
if do_AnalysisBasic == 1

    % --- --- --- --- --- --- --- --- --- --- ---
    % --- load results if they are not yet loaded
    % --- --- --- --- --- --- --- --- --- --- ---
    if ~exist('res')
        load([RESULTFOLDER  filesep iResultFile])
    end

    % disp('returning here XY'), return

    % --- first extract how many differnet distractor conditions and delayDurations were used:
    tmpDelays = [0.25:0.25:6.0];
    iConditions_nDistractors = [1 2 3];
    iConditions_delays = tmpDelays;
    iConditions = [];icnt=0;
    for j=1:length(iConditions_nDistractors)
        for k=1:length(iConditions_delays)
            icnt=icnt+1;
            iConditions(icnt,1:2) = [iConditions_nDistractors(j) iConditions_delays(k)];
        end
    end

    metrics_wm = [];
    metrics_wm.subject = [];
    metrics_wm.conditions_nDistr_delays = iConditions
    metrics_wm.allTrials = {};
    % --- averages
    metrics_wm.accuracyMeanSE = [];
    metrics_wm.RT_MeanSE = [];
    metrics_wm.accuracyMeanSE_LowTDS = [];
    metrics_wm.RT_MeanSE_LowTDS = [];
    metrics_wm.accuracyMeanSE_HighTDS = [];
    metrics_wm.RT_MeanSE_HighTDS = [];


    % --- --- --- --- --- --- --- --- ---
    % --- for each session:
    % --- --- --- --- --- --- --- --- ---
    for iD = 1:length(res)
if isempty(metrics_wm.subject) & ~isempty(res(iD).SubjectID)    
    metrics_wm.subject = res(iD).SubjectID{1};
end

        % --- --- --- --- --- --- --- --- ---
        % --- detemrine the delays used in this session:
        % --- --- --- --- --- --- --- --- ---
        res(iD).delayDurationOrig = res(iD).delayDuration;
        res(iD).delayDuration = round(res(iD).delayDuration .* 1000)/1000;
        for iTD=1:length(tmpDelays)
            sel = find( abs(res(iD).delayDuration - tmpDelays(iTD))<0.05 );
            res(iD).delayDuration(sel) = tmpDelays(iTD);
        end

        % --- --- --- --- --- --- --- --- ---
        % --- collect all trial outcomes for each combination of nDistr x Delay
        % --- --- --- --- --- --- --- --- ---
        for iM = 1:size(metrics_wm.conditions_nDistr_delays,1)
            sel = find(metrics_wm.conditions_nDistr_delays(iM,1) == res(iD).nDistractors ...
                & metrics_wm.conditions_nDistr_delays(iM,2) == res(iD).delayDuration )
            if isempty(sel), continue, end
            if length(metrics_wm.allTrials)<iM, metrics_wm.allTrials{iM} = []; end
            metrics_wm.allTrials{iM} = cat(1,metrics_wm.allTrials{iM}, [ res(iD).Outcome(sel)' ...
                res(iD).ReactionTime(sel)' ...
                repmat(iD,length(sel),1) ...
                res(iD).TDsharedFeatures(sel,1:metrics_wm.conditions_nDistr_delays(iM,1)) ]);
        end
    end

    % --- --- --- --- --- --- --- --- ---
    % --- compute condition meanAccuracy and meanRT
    % --- --- --- --- --- --- --- --- ---
    for iM = 1:size(metrics_wm.conditions_nDistr_delays,1)

        metrics_wm.accuracyMeanSE(iM,1:4) = NaN;
        metrics_wm.RT_MeanSE(iM,1:4) = NaN;
        metrics_wm.accuracyMeanSE_LowTDS(iM,1:4)  = NaN;
        metrics_wm.RT_MeanSE_LowTDS(iM,1:4)  = NaN;
        metrics_wm.accuracyMeanSE_HighTDS(iM,1:4)  = NaN;
        metrics_wm.RT_MeanSE_HighTDS(iM,1:4)  = NaN;
        if iM>length(metrics_wm.allTrials), continue, end
        if isempty(metrics_wm.allTrials{iM}), continue, end
        tmp = (metrics_wm.allTrials{iM}(:,1))>0;
        metrics_wm.accuracyMeanSE(iM,1:4) = [nanmean(tmp) nanmedian(tmp) nanstd(tmp)./sqrt(length(tmp)) length(tmp) ] ;
        tmpRT = metrics_wm.allTrials{iM}(:,2);
        metrics_wm.RT_MeanSE(iM,1:4) = [nanmean(tmpRT) nanmedian(tmpRT) nanstd(tmpRT)./sqrt(length(tmpRT)) length(tmpRT) ] ;

        tmpTDmean = nanmean(metrics_wm.allTrials{iM}(:,4:end),2)
        MMean = mean(nanmean(metrics_wm.allTrials{iM}(:,4:end)));

        selLow = find(tmpTDmean <= MMean);
        tmp = (metrics_wm.allTrials{iM}(selLow,1))>0; % >0 because token values are -1 and +2...
        metrics_wm.accuracyMeanSE_LowTDS(iM,1:4) = [nanmean(tmp) nanmedian(tmp) nanstd(tmp)./sqrt(length(tmp)) length(tmp) ] ;
        metrics_wm.RT_MeanSE_LowTDS(iM,1:4) = [nanmean(tmpRT(selLow)) nanmedian(tmpRT(selLow)) nanstd(tmpRT(selLow))./sqrt(length(tmpRT(selLow))) length(tmpRT(selLow)) ] ;

        selHigh = find(tmpTDmean > MMean);
        tmp = (metrics_wm.allTrials{iM}(selHigh,1))>0; % >0 because token values are -1 and +2...
        metrics_wm.accuracyMeanSE_HighTDS(iM,1:4) = [nanmean(tmp) nanmedian(tmp) nanstd(tmp)./sqrt(length(tmp)) length(tmp) ] ;
        metrics_wm.RT_MeanSE_HighTDS(iM,1:4) = [nanmean(tmpRT(selHigh)) nanmedian(tmpRT(selHigh)) nanstd(tmpRT(selHigh))./sqrt(length(tmpRT(selHigh))) length(tmpRT(selHigh)) ] ;

    end
    % disp('returning here XY'), return

    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    % --- save analysis/metrics results in MAT file
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    save([RESULTFOLDER filesep iResultFileMetrics],'res','metrics_wm','-V7.3')
    disp(sprintf('saved metrics mat-file of n=%d datasets in %s.',length(res), iResultFileMetrics)),

end




if do_AnalysisBasicPlot == 1

    % --- --- --- --- --- --- --- --- --- --- ---
    % --- load metrics results if they are not yet loaded
    % --- --- --- --- --- --- --- --- --- --- ---
    if ~exist('metrics_wm')
        load([RESULTFOLDER filesep iResultFileMetrics])
    end

    DO_SAVEFIGURES = 1;


   iSubject = metrics_wm.subject;
   nDistractorConditions = 3;
    nc=nDistractorConditions; nr=2; cs=1;
        figure('Color','w'), hold on, %        set(gcf,'Position',[1   800   579   201])
        set(gcf,'Position',[2    86   557   915])

   for iNDIS=1:nDistractorConditions
        sel = find(metrics_wm.conditions_nDistr_delays(:,1)==iNDIS);
        nDistractors = unique(metrics_wm.conditions_nDistr_delays(sel,1));
        iConditionLabel = sprintf('n%dDistr',nDistractors)
        iDelays = unique(metrics_wm.conditions_nDistr_delays(sel,2));
        for j=1:length(iDelays)
            tmp = find(metrics_wm.conditions_nDistr_delays(sel,2)==iDelays(j));
            if isempty(tmp), continue, end
            iM = sel(tmp);
            iZ1(j,:) = metrics_wm.accuracyMeanSE(iM,1:4);
            iZ2(j,:) = metrics_wm.accuracyMeanSE_LowTDS(iM,1:4);
            iZ3(j,:) = metrics_wm.accuracyMeanSE_HighTDS(iM,1:4);
            iR1(j,:) = metrics_wm.RT_MeanSE(iM,1:4);
            iR2(j,:) = metrics_wm.RT_MeanSE_LowTDS(iM,1:4);
            iR3(j,:) = metrics_wm.RT_MeanSE_HighTDS(iM,1:4);
        end

        t = find(~isnan(iZ1(:,1)));
        if isempty(t), continue, end

        iX = iDelays(t);
        iMin = 0.5;
        %  --- accuracy
        subplot(nc,nr,cs), hold on, cs=cs+1;
        iZ = iZ1(t,1);    iZSE = iZ1(t,3); iMin = min(([min(iZ-iZSE) iMin]));
        errorbar(iX, iZ,iZSE,'o-','Color','k','MarkerFaceColor','k'), hold on

        iZ = iZ2(t,1);    iZSE = iZ2(t,3); iMin = min(([min(iZ-iZSE) iMin]));
        errorbar(iX, iZ,iZSE,'o--','Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]), hold on
        iZ = iZ3(t,1);    iZSE = iZ3(t,3); iMin = min(([min(iZ-iZSE) iMin]));
        errorbar(iX, iZ,iZSE,'o-.','Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]), hold on

        set(gca,'xlim',[0 max(iX)+0.25],'tickdir','out','xtick',[0:0.25:max(iX)+0.25])
        ylabel('Accuracy'), xlabel('delay (s)'),
        title(sprintf('%s.min%d.max%d.%s',iSubject, min(iZ1(t,4)), max(iZ1(t,4)),iConditionLabel))
        set(gca,'ylim',[iMin 1])
        plot(get(gca,'xlim'),[0.75 0.75],'Color',[0.7 0.7 0.7])
        plot(get(gca,'xlim'),[0.5 0.5],'Color',[0.7 0.7 0.7])
        plot(get(gca,'xlim'),[0.25 0.25],'Color',[0.7 0.7 0.7])
        legend('avg','lowSim','highSim','location','Southwest')

        %  --- rt
        subplot(nc,nr,cs), hold on, cs=cs+1;
        iR = iR1(t,1);    iRSE = iR1(t,3);
        errorbar(iX, iR,iRSE,'o-','Color','k','MarkerFaceColor','k'), hold on

        iR = iR2(t,1);    iRSE = iR2(t,3);
        errorbar(iX, iR,iRSE,'o--','Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]), hold on
        iR = iR3(t,1);    iRSE = iR3(t,3);
        errorbar(iX, iR,iRSE,'o-.','Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6]), hold on

        set(gca,'xlim',[0 max(iX)+0.25],'tickdir','out','xtick',[0:0.25:max(iX)+0.25])
        ylabel('RespT'), xlabel('delay (s)'),
        title(sprintf('%s.min%d.max%d',iSubject, min(iR1(t,4)), max(iR1(t,4)),iConditionLabel))
        %set(gca,'ylim',[0.3 1])
        %plot(get(gca,'xlim'),[0.75 0.75],'Color',[0.7 0.7 0.7])
        %plot(get(gca,'xlim'),[0.5 0.5],'Color',[0.7 0.7 0.7])
        %plot(get(gca,'xlim'),[0.25 0.25],'Color',[0.7 0.7 0.7])
        legend('avg','lowSim','highSim','location','Southwest')

    end

        if DO_SAVEFIGURES
            currentDateNumber = now;
            currentDateNumber = datestr(currentDateNumber, 'yyyymmdd');
            figurefilename = sprintf('fig_wm_AcrossSessions_type1_%s_%s.pdf',iMonkey,currentDateNumber);
            saveas(gcf, [FIGURE_Folder filesep figurefilename], 'pdf');
        end
        close gcf


end


disp('end of analyze_visualSearch.m'), return
