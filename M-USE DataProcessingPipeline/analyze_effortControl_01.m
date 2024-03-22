

%
% This script contains some basic analyzes of the Effort Control task.
% It assumes that the script analyze_preprocessSessions_01.m was run to
% save trialdata from 1 or more sessions containing the EC task

%
% Copyright (c) 2023, Marcus Watson, Thilo Womelsdorf
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
%

% id 2023-08-09 analyze_effortControl_01
% --- --- --- --- --- --- --- ---
% --- folder with preprocecced data
% --- --- --- --- --- --- --- ---

do_LoadData = 1;
do_AnalysisBasic = 1;
do_AnalysisBasicPlot = 1;

% --- --- --- --- --- --- --- ---
% --- folder with preprocecced data
% --- --- --- --- --- --- --- ---
HOME_FOLDER = ['/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS' filesep];
%'/Users/thiwom/Main/projects/201912_P_KioskFluBehavior/M_USE/GAMES_MATLAB_ANALYSIS'
MUSEMATFOLDERNames = {}; iResultFolder = ''; iResultFile = '';


% SessionID_WM = {'Wotan_EC_all_01'}

%SessionID_EC = {'Frey_EC_all_01'}
SessionID_EC = {'Wotan_EC_all_01'}

% MUSEMATFOLDERNames = 'MUSEMAT01_VS_EC_FL_MZG_Frey'
if strcmp(SessionID_EC{1},'Frey_EC_all_01')
    MUSEMATFOLDERNames{1} = 'MUSEMAT01_WM_EC_FL_MZG_Frey'
    % MUSEMATFOLDERNames{2} = 'MUSEMAT01_VS_EC_FL_MRT_Frey'
    iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
    iResultFile   = 'EC_Frey_01_20240201'; %'WM_01_20231205';%'WM_01_20230907';
end
if strcmp(SessionID_EC{1},'Wotan_EC_all_01')
    MUSEMATFOLDERNames{1} =  'MUSEMAT01_WM_EC_FL_MZG_Wotan';
    %MUSEMATFOLDERNames{2} =  'MUSEMAT01_VS_EC_FL_MRT_Wotan';
    iResultFolder = [MUSEMATFOLDERNames{1} '_RES'];
    iResultFile   = 'EC_Wotan_01_20240201';%'WM_01_20230907';
end

iMonkey = SessionID_EC{1}(1:findstr(SessionID_EC{1},'_')-1);

% --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
% --- The files generates by this script and the folder name for them:
iResultFileMetrics = [iResultFile '_Metrics'];
RESULTFOLDER = [HOME_FOLDER filesep iResultFolder '_MAT']; if ~exist(RESULTFOLDER), mkdir(RESULTFOLDER),end
METRICSFOLDER  = [RESULTFOLDER '_METRICS']; if ~exist(METRICSFOLDER), mkdir(METRICSFOLDER),end
FIGURE_Folder = [RESULTFOLDER '_FIG'];    if ~exist(FIGURE_Folder), mkdir(FIGURE_Folder), end

if do_LoadData==1
    % iResultFile = ['DAT01_' Monkey '_' EXP_ID '_' iSessionName];

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
    res.nTrials             = [];
    res.effort_token_cnd    = [];
    res.choice_nTouches     = [];
    res.choice_side         = [];
    res.choice_effort_reward = [];
    res.choiceDuration      = [];
    % trialData(iT).LoadTrialStims_EndTimeRelative - trialData(iT).LoadTrialStims_StartTimeRelative;
    res.heldTooShortLong    = [];

    % --- for analysis of overall response options
    res.Eff_X = [];
    res.Rew_X = [];
    res.ChosenEffort_prop_nn = [];
    res.ChosenReward_prop_nn = [];

    % --- for analysis of overall choices relative to the difference in
    % --- effort/reward
    res.effort_token_diff = [];
    res.EffDiff_X = [];
    res.RewDiff_X = [];
    res.ChosenHigherRew_EffDiff__prop_nn = [];
    res.ChosenHigherEff_RewDiff__prop_nn = [];
    res.ChosenHigherCMB_EffRew_prophigherEff_nn = [];
    res.ChosenHigherCMB_EffRew_prophigherRew_nn = [];

    for iD=1:length(datasets)

        res(iD).dataset = [datasets{iD}];
        res(iD).datafolder = [datasetsDIR{iD}];
        res(iD).log = [];

        % DAT01_Frey_VS_EC_FL_MRT_Session_08_07_23__12_39_49__SubjectID_Frey
        disp(sprintf('processing iD=%d, %s',iD,datasets{iD})),
        %in = load([DATAFOLDER filesep datasets{iD}]); disp('.')
        in = load([datasetsDIR{iD} filesep datasets{iD}]); disp('.')

        % idx_all = find(contains(in.dat.taskLabel,'FlexLearning'));
        % idx_all = find(contains(in.dat.taskLabel,'VisualSearch'));
        idx_all = find(contains(in.dat.taskLabel,'EffortControl'));
        % idx_all = find(contains(in.dat.taskLabel,'MazeGuided'));
        % idx_all = find(contains(in.dat.taskLabel,'Maze'));
        %idx_all = find(contains(in.dat.taskLabel,'WorkingMemory'));
        %idx = strmatch('EffortControl',in.dat.taskLabel);
        if isempty(idx_all), res(iD).log=[]; res(iD).log{end+1} = sprintf('task isempty...iD=%d, %s',iD,datasets{iD}); continue, end
        
        for iIDX = 1:length(idx_all)
            idx = idx_all(iIDX);

            if isempty(idx), continue, end
            if idx > length(in.dat.trialData), res(iD).log=[]; res(iD).log{end+1} = sprintf('no trials no task...iD=%d, %s',iD,datasets{iD}); continue, end

            trialData = in.dat.trialData{idx};

            res(iD).nTrials             = length(trialData);

            for iT=1:res(iD).nTrials
                res(iD).blockNum(iT) = trialData(iT).BlockCount;

                res(iD).effort_token_cnd(iT,:) = [  trialData(iT).ClicksNeededLeft,...
                    trialData(iT).ClicksNeededRight,...
                    trialData(iT).NumCoinsLeft,...
                    trialData(iT).NumCoinsRight ];

                if ~isfield(trialData, 'ClicksPerOutline') & isfield(trialData,'ClicksNeeded')
                    trialData(iT).ClicksPerOutline = trialData(iT).ClicksNeeded;
                end

                res(iD).choice_nTouches(iT) = trialData(iT).ClicksPerOutline;
                res(iD).choice_side(iT)     = strcmp(trialData(iT).ChosenSide,'Right'); % Left is 0, Right is 1
                chosenEffort = NaN;
                if strcmp(trialData(iT).ChosenEffort,'Same'), chosenEffort=0;
                elseif strcmp(trialData(iT).ChosenEffort,'Lower'), chosenEffort=1;
                elseif strcmp(trialData(iT).ChosenEffort,'Higher'), chosenEffort=2;
                end
                chosenReward = NaN;
                if strcmp(trialData(iT).ChosenReward,'Same'), chosenReward=0;
                elseif strcmp(trialData(iT).ChosenReward,'Lower'), chosenReward=1;
                elseif strcmp(trialData(iT).ChosenReward,'Higher'), chosenReward=2;
                end
                res(iD).choice_effort_reward_label = '0:Same,1:Lower,2:Higher';
                res(iD).choice_effort_reward(iT,:) = [chosenEffort chosenReward];
                res(iD).choiceDuration(iT)  = trialData(iT).TimeTakenToChoose;
                % trialData(iT).LoadTrialStims_EndTimeRelative - trialData(iT).LoadTrialStims_StartTimeRelative;

                res(iD).heldTooShortLong(iT,:)  = [trialData(iT).HeldTooShort trialData(iT).HeldTooLong];

                sprintf('done with %d trials',iT)
            end
        end
    end
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    % --- save analysis results in RES file
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    save([RESULTFOLDER filesep iResultFile],'res','-V7.3')
    fprintf('saved analysis of n=%d datasets in %s \n',iD, iResultFile),

end


if do_AnalysisBasic == 1

    % --- --- --- --- --- --- --- --- --- --- ---
    % --- load results if they are not yet loaded
    % --- --- --- --- --- --- --- --- --- --- ---
    if ~exist('res')
        load([RESULTFOLDER iResultFile])
    end


    metrics_ec = [];
    for iD = 1:length(res)
        metrics_ec(iD).dataset = [];
        metrics_ec(iD).Eff_X = [];
        metrics_ec(iD).Rew_X = [];
        metrics_ec(iD).ChosenEffort_prop_nn = [];
        metrics_ec(iD).ChosenReward_prop_nn = [];
        metrics_ec(iD).effort_token_diff = [];
        metrics_ec(iD).EffDiff_X = [];
        metrics_ec(iD).RewDiff_X = [];
        metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn = [];
        metrics_ec(iD).ChosenHigherEff_RewDiff__prop_nn = [];
        metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn = [];
        metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn = [];
        metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn = [];
        metrics_ec(iD).ChosenHigherEff_RewDiff__prop_nn = [];
        metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn = [];
        metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn = [];
        metrics_ec(iD).XDiff_EFF_REW = [];

        % disp('rrr'), return
        if isempty(res(iD).effort_token_cnd), continue, end

        % --- sort into the specific Reward and Effort conditions:
        metrics_ec(iD).Eff_X = unique(res(iD).effort_token_cnd(:,1:2));
        metrics_ec(iD).Rew_X = unique(res(iD).effort_token_cnd(:,3:4));

        metrics_ec(iD).ChosenEffort_prop_nn = [];
        metrics_ec(iD).ChosenReward_prop_nn = [];

        metrics_ec(iD).dataset =  res(iD).dataset;
        metrics_ec(iD).effort_token_cnd = res(iD).effort_token_cnd;



        % --- how often was each effort level chosen ?
        for k=1:length(metrics_ec(iD).Eff_X)
            % --- how often was this effort level shown
            idx_n1   = find( res(iD).effort_token_cnd(:,1) == metrics_ec(iD).Eff_X(k));
            idx_n2   = find( res(iD).effort_token_cnd(:,2) == metrics_ec(iD).Eff_X(k));

            % --- remove trials when both were the same
            [a,b,c]=intersect(idx_n1,idx_n2);
            if ~isempty(b), idx_n1(b)=[]; idx_n2(c)=[]; end

            % --- how often was this effort level chosen when it was shown
            n1chosen = find( res(iD).choice_side(idx_n1) == 0);  % chosen left
            n2chosen = find( res(iD).choice_side(idx_n2) == 1);  % chosen right

            nShown = length(idx_n1) + length(idx_n2);
            nChosen = length(n1chosen) + length(n2chosen);

            metrics_ec(iD).ChosenEffort_prop_nn(k,:) = [(nChosen / nShown)  nChosen  nShown];
        end

        % --- how often was each reward level chosen ?
        for k=1:length(metrics_ec(iD).Rew_X)
            % --- how often was this effort level shown
            idx_n1   = find( res(iD).effort_token_cnd(:,3) == metrics_ec(iD).Rew_X(k));
            idx_n2   = find( res(iD).effort_token_cnd(:,4) == metrics_ec(iD).Rew_X(k));

            % --- remove trials when both were the same
            [a,b,c]=intersect(idx_n1,idx_n2);
            if ~isempty(b), idx_n1(b)=[]; idx_n2(c)=[]; end

            % --- how often was this effort level chosen when it was shown
            n1chosen = find( res(iD).choice_side(idx_n1) == 0);  % chosen left
            n2chosen = find( res(iD).choice_side(idx_n2) == 1);  % chosen right

            nShown = length(idx_n1) + length(idx_n2);
            nChosen = length(n1chosen) + length(n2chosen);

            metrics_ec(iD).ChosenReward_prop_nn(k,:) = [(nChosen / nShown)  nChosen  nShown];
        end


        % --- --- --- --- --- --- --- --- --- ---
        % --- how likely is higher reward chosen given each effort-difference ?
        % --- how likely is higher effort level chosen given each reward-difference ?
        % --- --- --- --- --- --- --- --- --- ---
        metrics_ec(iD).effort_token_diff = [abs(diff(res(iD).effort_token_cnd(:,1:2)')') abs(diff(res(iD).effort_token_cnd(:,3:4)')')];

        % [conditionsDiff,b,c] = unique(res(iD).effort_token_diff(:,1:2),'rows')

        metrics_ec(iD).EffDiff_X = unique(metrics_ec(iD).effort_token_diff(:,1));
        metrics_ec(iD).RewDiff_X = unique(metrics_ec(iD).effort_token_diff(:,2));

        metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn = nan(length(metrics_ec(iD).EffDiff_X),3);
        metrics_ec(iD).ChosenHigherEff_RewDiff__prop_nn = nan(length(metrics_ec(iD).RewDiff_X),3);

        metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn = nan(length(metrics_ec(iD).EffDiff_X) * length(metrics_ec(iD).RewDiff_X), 3);
        metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn = nan(length(metrics_ec(iD).EffDiff_X) * length(metrics_ec(iD).RewDiff_X), 3);

        metrics_ec(iD).XDiff_EFF_REW = nan(length(metrics_ec(iD).EffDiff_X) * length(metrics_ec(iD).RewDiff_X));

        icnt = 0;
        for k=1:length(metrics_ec(iD).EffDiff_X)
            for j=1:length(metrics_ec(iD).RewDiff_X)

                if j==1
                    % --- how often was this effort-level difference shown
                    idxk = find(metrics_ec(iD).effort_token_diff (:,1) == metrics_ec(iD).EffDiff_X(k));

                    % --- how often was higher reward chosen
                    nHigherorLower = length(find(res(iD).choice_effort_reward(idxk,2) == 1 | res(iD).choice_effort_reward(idxk,2) == 2));
                    nHigher =length( find(res(iD).choice_effort_reward(idxk,2) == 2));
                    metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn(k,:) = [nHigher/nHigherorLower nHigher nHigherorLower];
                end
                if k==1
                    % --- how often was this reward level difference shown
                    idxj = find(metrics_ec(iD).effort_token_diff(:,2) == metrics_ec(iD).RewDiff_X(j));

                    % --- how often was higher effort chosen
                    nHigherorLower = length( find(res(iD).choice_effort_reward(idxj,1) == 1 | res(iD).choice_effort_reward(idxj,1) == 2));
                    nHigher = length( find(res(iD).choice_effort_reward(idxj,1) == 2));
                    metrics_ec(iD).ChosenHigherEff_RewDiff__prop_nn(j,:) = [nHigher/nHigherorLower  nHigher nHigherorLower];
                end

                % how often did this combination occur?
                idx = find( (metrics_ec(iD).effort_token_diff(:,1) == metrics_ec(iD).EffDiff_X(k))   &   (metrics_ec(iD).effort_token_diff(:,2) == metrics_ec(iD).RewDiff_X(j))  );

                % how often was the higher effort /higher reward  chosen with this combination of differences ?
                % nHigherorLower = length(find(res(iD).choice_effort_reward(idxk,2) == 1 | res(iD).choice_effort_reward(idxk,2) == 2));
                nHigherEff =length( find(res(iD).choice_effort_reward(idx,1) == 2));
                nHigherRew =length( find(res(iD).choice_effort_reward(idx,2) == 2));

                icnt = icnt + 1;
                metrics_ec(iD).X_EFF_REW(icnt,1:2) = [ metrics_ec(iD).EffDiff_X(k) metrics_ec(iD).RewDiff_X(j) ];
                metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn(icnt,1:3) = [nHigherEff/length(idx) nHigherEff length(idx) ];
                metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn(icnt,1:3) = [nHigherRew/length(idx) nHigherRew length(idx) ];
            end
        end
    end

    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    % --- save analysis/metrics results in MAT file
    % --- --- --- --- --- --- --- --- --- --- --- --- ---
    save([RESULTFOLDER filesep iResultFileMetrics],'res','metrics_ec','-V7.3')
    disp(sprintf('saved metrics mat-file of n=%d datasets in %s.',length(res), iResultFileMetrics)),

end % end of analysis


if do_AnalysisBasicPlot == 1


    % --- --- --- --- --- --- --- --- --- --- ---
    % --- load metrics results if they are not yet loaded
    % --- --- --- --- --- --- --- --- --- --- ---
    DO_SAVEFIGURES = 1;
    PLOT_INDIVIDUALSESSIONS = 1;
    PLOT_ACROSSSESSIONS = 1;

    if ~exist('metrics_ec')
        load([RESULTFOLDER filesep iResultFileMetrics])
    end
    % --- --- --- --- --- --- --- --- --- --- ---
    % --- Plot individual sessions
    % --- --- --- --- --- --- --- --- --- --- ---
    if PLOT_INDIVIDUALSESSIONS==1
        for iD = 1: length(metrics_ec)
            %iD = 15; % session id

            iSessionName = metrics_ec(iD).dataset;
            if isempty(metrics_ec(iD).Eff_X), continue, end
            MS=5;
            figure('Color','w'), hold on,
            % set(gcf,'Position',[3   388   560   613])%4   581   560   420
            set(gcf,'Position',[2    86   557   915])


            nr=4;nc=2;cnt=1;
            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            textboxTxT = sprintf('%s\n%s',iSessionName(1:findstr(iSessionName,'Session')-2), ...
                iSessionName(findstr(iSessionName,'Session'):end));
            % Create an annotation textbox
            annotation('textbox', [0.1, 0.1, 0.8, 0.8], 'String', textboxTxT, 'FitBoxToText', 'on', 'EdgeColor', 'none');
            title('effort control session');
            axis off; % Turn off the axis if you don't need it

            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            isel = find(~isnan(metrics_ec(iD).Eff_X));
            x = metrics_ec(iD).Eff_X(isel);
            y = metrics_ec(iD).ChosenEffort_prop_nn(isel,1);
            [out] = get_fitValues_EC(x,y)
            if any(y>0)
                plot(out.x, out.y, 'ro', 'MarkerSize',MS); hold on,
                plot(out.xFit, out.yFit, '-r'), hold on,
                set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
                set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
                plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
                plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
            end
            xlabel('numEffort'), ylabel('prop chosen'), title(sprintf('y50 %.2f,x50 %.2f p chose higher Eff',out.yValueAtSlopePeak,out.fitParams(1)))
            axis square

            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            %plot(metrics_ec(iD).Rew_X, metrics_ec(iD).ChosenReward_prop_nn(:,1),'bo-', 'MarkerSize',MS);
            isel = find(~isnan(metrics_ec(iD).Rew_X));
            x = metrics_ec(iD).Rew_X(isel);
            y = metrics_ec(iD).ChosenReward_prop_nn(isel,1),
            [out] = get_fitValues_EC(x,y,'positive')
            if any(y>0)
                plot(out.x, out.y, 'ro', 'MarkerSize',MS); hold on,
                plot(out.xFit, out.yFit, '-r'), hold on,
                set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
                set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
                plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
                plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
            end
            xlabel('numTokens'), ylabel('prop chosen'), title(sprintf('y50 %.2f,x50 %.2f',out.yValueAtSlopePeak,out.fitParams(1)))
            axis square

            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            %plot(metrics_ec(iD).EffDiff_X, metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn(:,1),'ro-', 'MarkerSize',MS);
            isel = find(~isnan(metrics_ec(iD).EffDiff_X));
            x = metrics_ec(iD).EffDiff_X(isel);
            y = metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn(isel,1);
            [out] = get_fitValues_EC(x,y)
            if any(y>0)
                plot(out.x, out.y, 'ro', 'MarkerSize',MS); hold on,
                plot(out.xFit, out.yFit, '-r'), hold on,
                set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
                set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
                plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
                plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
            end
            xlabel('Eff Diff'), ylabel('p chosing higher rew'), title(sprintf('y50 %.2f,x50 %.2f',out.yValueAtSlopePeak,out.fitParams(1)))
            axis square

            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            isel = find(~isnan(metrics_ec(iD).RewDiff_X));
            x=metrics_ec(iD).RewDiff_X(isel);
            y=metrics_ec(iD).ChosenHigherEff_RewDiff__prop_nn(isel,1);
            %xlabel('rew Diff'), ylabel('p chose higher effort')
            [out] = get_fitValues_EC(x,y)
            plot(out.x, out.y, 'ro', 'MarkerSize',MS); hold on,
            plot(out.xFit, out.yFit, '-b'), hold on,
            if any(y>0)
                set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
                set(gca,'ylim', [min([min(out.yFit) min(y)]) .* 0.8  max(y).*1.1] )
                plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
                plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
            end
            xlabel('Rew Diff'), ylabel('p chose higher effort'), title(sprintf('y50 %.2f,x50 %.2f',out.yValueAtSlopePeak,out.fitParams(1)))
            axis square

            % nr=3;nc=2;cnt=1;
            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            XYZ = [ metrics_ec(iD).X_EFF_REW metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn(:,3) ];
            [iX, iY , iZ] = get_X_Y_Zmatrix(XYZ);
            imagesc(iX,iY,iZ), axis tight, colorbar,
            set(gca, 'YDir', 'reverse')
            ylabel('diffEff'), xlabel('diffRew'), title('n cmb')
            axis square


            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            XYZ = [ metrics_ec(iD).X_EFF_REW metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn(:,1) ];
            [iX, iY, iZ] = get_X_Y_Zmatrix(XYZ);
            imagesc(iY,iX,iZ), axis tight, colorbar,
            set(gca, 'YDir', 'reverse')
            ylabel('diffEff'), xlabel('diffRew'), title('p chose higher Eff')
            axis square

            subplot(nr,nc,cnt), hold on,cnt=cnt+1;
            XYZ = [ metrics_ec(iD).X_EFF_REW metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn(:,1) ];
            [iX ,iY, iZ] = get_X_Y_Zmatrix(XYZ);
            imagesc(iY,iX,iZ), axis tight, colorbar,
            set(gca, 'YDir', 'reverse')
            ylabel('diffEff'), xlabel('diffRew'), , title('p chose higher Rew')
            axis square


            if DO_SAVEFIGURES
                currentDateNumber = now;
                currentDateNumber = datestr(currentDateNumber, 'yyyymmdd');
                figurefilename = sprintf('fig_ec_session%2.d_type1_%s_%s.pdf',iD,iMonkey, currentDateNumber);
                saveas(gcf, [FIGURE_Folder filesep figurefilename], 'pdf');
            end
            close gcf
        end
    end

    %disp(' XX adjust metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn'), return



    if PLOT_ACROSSSESSIONS==1

        metrics_ec_sum = [];
        N = length(metrics_ec);
        metrics_ec_sum.EFF_X=nan(N,4);
        metrics_ec_sum.EFF_Y=nan(N,4);
        metrics_ec_sum.EFF_fitParams=nan(N,4);
        metrics_ec_sum.EFF_yValueAtSlopePeak=nan(N,1);
        metrics_ec_sum.EFFDiff_X=nan(N,4);
        metrics_ec_sum.EFFDiff_Y=nan(N,4);
        metrics_ec_sum.EFFDiff_fitParams=nan(N,4);
        metrics_ec_sum.EFFDiff_yValueAtSlopePeak=nan(N,1);

        metrics_ec_sum.REW_X=nan(N,4);
        metrics_ec_sum.REW_Y=nan(N,4);
        metrics_ec_sum.REW_fitParams=nan(N,4);
        metrics_ec_sum.REW_yValueAtSlopePeak=nan(N,1);

        metrics_ec_sum.REWDiff_X=nan(N,4);
        metrics_ec_sum.REWDiff_Y=nan(N,4);
        metrics_ec_sum.REWDiff_fitParams=nan(N,4);
        metrics_ec_sum.REWDiff_yValueAtSlopePeak=nan(N,1);

        metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn = [];
        metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherRew_nn = [];

        % first average accross sessions:
        iS=0;
        for iD = 1: N
            if isempty(metrics_ec(iD).Eff_X), continue, end
            sel = find(~isnan(metrics_ec(iD).Eff_X));
            if length(sel)<=2, continue, end

            iSessionName = metrics_ec(iD).dataset;

            iS=iS+1;

            x = metrics_ec(iD).Eff_X(sel);
            y = metrics_ec(iD).ChosenEffort_prop_nn(sel,1);
            [out] = get_fitValues_EC(x,y);
            metrics_ec_sum.EFF_X(iS,1:length(x)) =x;
            metrics_ec_sum.EFF_Y(iS,1:length(y)) =y;
            metrics_ec_sum.EFF_fitParams(iS,1:length(out.fitParams)) = out.fitParams;
            metrics_ec_sum.EFF_yValueAtSlopePeak(iS) =out.yValueAtSlopePeak;
            % sprintf('y50 %.2f,x50 %.2f p chose higher Eff',out.yValueAtSlopePeak,out.fitParams(1))

            x = metrics_ec(iD).Rew_X(sel);
            y = metrics_ec(iD).ChosenReward_prop_nn(sel,1);
            [out] = get_fitValues_EC(x,y);
            metrics_ec_sum.REW_X(iS,1:length(x)) =x;
            metrics_ec_sum.REW_Y(iS,1:length(y)) =y;
            metrics_ec_sum.REW_fitParams(iS,1:length(out.fitParams)) = out.fitParams;
            metrics_ec_sum.REW_yValueAtSlopePeak(iS) =out.yValueAtSlopePeak;
            % sprintf('y50 %.2f,x50 %.2f p chose higher Eff',out.yValueAtSlopePeak,out.fitParams(1))

            sel = find(~isnan(metrics_ec(iD).EffDiff_X));
            if length(sel)>2,
                x = metrics_ec(iD).EffDiff_X(sel);;
                y = metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn(sel,1);
                [out] = get_fitValues_EC(x,y);
                metrics_ec_sum.EFFDiff_X(iS,1:length(x))  = x;
                metrics_ec_sum.EFFDiff_Y(iS,1:length(y))  = y;
                metrics_ec_sum.EFFDiff_fitParams(iS,1:length(out.fitParams)) = out.fitParams;
                metrics_ec_sum.EFFDiff_yValueAtSlopePeak(iS) =out.yValueAtSlopePeak;
            end

            sel = find(~isnan(metrics_ec(iD).RewDiff_X));
            if length(sel)>2,
                x = metrics_ec(iD).RewDiff_X(sel);
                y = metrics_ec(iD).ChosenHigherEff_RewDiff__prop_nn(sel,1);
                [out] = get_fitValues_EC(x,y);
                metrics_ec_sum.REWDiff_X(iS,1:length(x)) = x;
                metrics_ec_sum.REWDiff_Y(iS,1:length(y))  = y;
                metrics_ec_sum.REWDiff_fitParams(iS,1:length(out.fitParams)) = out.fitParams;
                metrics_ec_sum.REWDiff_yValueAtSlopePeak(iS) =out.yValueAtSlopePeak;

                %metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn = nan(N,1);
                %metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherRew_nn = nan(N,1);
                %ChosenHigherCMB_EffRew_prophigherEff_nn
                %ChosenHigherCMB_EffRew_prophigherRew_nn
            end

            XYZ = [ metrics_ec(iD).X_EFF_REW metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn(:,1:3) ];
            metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn = cat(1, metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn, XYZ);
            XYZ = [ metrics_ec(iD).X_EFF_REW metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherRew_nn(:,1:3) ];
            metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherRew_nn = cat(1, metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherRew_nn, XYZ);

        end

        %metrics_ec(iD).ChosenHigherCMB_EffRew_prophigherEff_nn
        iMatrix = metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn;
        dsel = [isnan(iMatrix(:,1)) | isnan(iMatrix(:,2))]; iMatrix(dsel,:)=[];
        [uniqueRows, ~, ic] = unique(iMatrix(:, 1:2), 'rows');
        sumChoice = accumarray(ic, iMatrix(:, end-1));
        sumTotal = accumarray(ic, iMatrix(:, end));
        XYZeffN = [uniqueRows sumTotal];
        XYZeff = [uniqueRows (sumChoice ./ sumTotal)];

        iMatrix = metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherRew_nn;
        dsel = [isnan(iMatrix(:,1)) | isnan(iMatrix(:,2))]; iMatrix(dsel,:)=[];
        [uniqueRows, ~, ic] = unique(iMatrix(:, 1:2), 'rows');
        sumChoice = accumarray(ic, iMatrix(:, end-1));
        sumTotal = accumarray(ic, iMatrix(:, end));
        XYZrewN = [uniqueRows sumTotal];
        XYZrew = [uniqueRows (sumChoice ./ sumTotal)];


        MS=5;
        figure('Color','w'), hold on,
            set(gcf,'Position',[2    86   557   915])

        nr=4;nc=2;cnt=1;
        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        %XYZ = [ metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn(:,[1 2 5]) ];
        [iX iY, iZ] = get_X_Y_Zmatrix(XYZeffN);
        imagesc(iY,iX,iZ), axis tight, colorbar,
        set(gca, 'YDir', 'reverse')
        ylabel('diffEff'), xlabel('diffRew'), title(sprintf('n cmb %s',iMonkey))
        axis square

        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        %XYZ = [  metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherEff_nn(:,[1 2 3]) ];
        [iX iY, iZ] = get_X_Y_Zmatrix(XYZeff);
        imagesc(iY,iX,iZ), axis tight, colorbar,
        set(gca, 'YDir', 'reverse')
        ylabel('diffEff'), xlabel('diffRew'), title('p chose higher Eff')
        axis square

        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        %XYZ = [ metrics_ec_sum.ChosenHigherCMB_EffRew_prophigherRew_nn(:,[1 2 3]) ];
        [iX iY, iZ] = get_X_Y_Zmatrix(XYZrew);
        %imagesc(iZ), axis tight, colorbar,
        imagesc(iY,iX,iZ), axis tight, colorbar,
        set(gca, 'YDir', 'reverse')
        ylabel('diffEff'), xlabel('diffRew'), , title('p chose higher Rew')
        axis square, 
  



        %subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        %nSessions = size(metrics_ec_sum.EFF_fitParams,1);
        %textboxTxT = sprintf('nSessions: %d\n',nSessions);
        %% Create an annotation textbox
        %annotation('textbox', [0.1, 0.1, 0.8, 0.8], 'String', textboxTxT, 'FitBoxToText', 'on', 'EdgeColor', 'none');
        %title('effort control across sessions');
        %axis off; % Turn off the axis if you don't need it


        %cnt = cnt+1; % skip one

        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        N = size(metrics_ec_sum.EFF_fitParams,1);
        x = metrics_ec_sum.EFF_X(1,:);
        y = nanmean(metrics_ec_sum.EFF_Y,1);
        yse = nanstd(metrics_ec_sum.EFF_Y)./N;
        [out] = get_fitValues_EC(x,y)
        errorbar(x, y, yse, 'ro-', 'MarkerSize',MS); hold on
        plot(out.xFit, out.yFit, '-r'), hold on,
        set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
        set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
        plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
        plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
        xlabel('numEffort'), ylabel('prop E chosen'), title(sprintf('y50 %.2f,x50 %.2f p chose higher Eff',out.yValueAtSlopePeak,out.fitParams(1)))
        axis square

        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        N = size(metrics_ec_sum.REW_fitParams,1);
        x = metrics_ec_sum.REW_X(1,:);
        y = nanmean(metrics_ec_sum.REW_Y,1);
        yse = nanstd(metrics_ec_sum.REW_Y)./N;
        [out] = get_fitValues_EC(x,y,'positive')
        errorbar(x, y, yse, 'ro', 'MarkerSize',MS);
        plot(out.xFit, out.yFit, '-r'), hold on,
        set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
        set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
        plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
        plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
        xlabel('numTokens'), ylabel('prop R chosen'), title(sprintf('y50 %.2f,x50 %.2f',out.yValueAtSlopePeak,out.fitParams(1)))
        axis square

        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        %plot(metrics_ec(iD).EffDiff_X, metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn(:,1),'ro-', 'MarkerSize',MS);
        %x = metrics_ec(iD).EffDiff_X(1,:);
        %y = metrics_ec(iD).ChosenHigherRew_EffDiff__prop_nn(:,1);
        N = size(metrics_ec_sum.EFFDiff_X,1);
        x = metrics_ec_sum.EFFDiff_X(1,:);
        y = nanmean(metrics_ec_sum.EFFDiff_Y,1);
        yse = nanstd(metrics_ec_sum.EFFDiff_Y)./N;
        [out] = get_fitValues_EC(x,y)
        plot(out.x, out.y, 'ro', 'MarkerSize',MS); hold on,
        plot(out.xFit, out.yFit, '-r'), hold on,
        set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
        set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
        plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
        plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
        xlabel('Eff Diff'), ylabel('p chosing higher rew'), title(sprintf('y50 %.2f,x50 %.2f',out.yValueAtSlopePeak,out.fitParams(1)))
        axis square

        %metrics_ec_sum.EFFDiff_fitParams(iS,1:length(out.fitParams)) = out.fitParams;
        %metrics_ec_sum.EFFDiff_yValueAtSlopePeak(iS) =out.yValueAtSlopePeak;

        subplot(nr,nc,cnt), hold on,cnt=cnt+1;
        N = size(metrics_ec_sum.REWDiff_X,1);
        x = metrics_ec_sum.REWDiff_X(1,:);
        y = nanmean(metrics_ec_sum.REWDiff_Y,1);
        yse = nanstd(metrics_ec_sum.REWDiff_Y)./N,

        %xlabel('rew Diff'), ylabel('p chose higher effort')
        [out] = get_fitValues_EC(x,y)
        errorbar(x, y, yse, 'ro', 'MarkerSize',MS); hold on,
        plot(out.xFit, out.yFit, '-b'), hold on,
        set(gca,'tickdir','out','xlim', [min(x)-2 max(x)+2 ],'xtick',[x])
        set(gca,'ylim', [min([min(out.yFit) min(y)]).*0.8  max(y).*1.1] )
        plot([min(get(gca,'xlim')) out.fitParams(1)], [repmat(out.yValueAtSlopePeak,1,2)], '-k'), hold on,
        plot([out.fitParams(1) out.fitParams(1)], [min(get(gca,'ylim')) out.yValueAtSlopePeak], '-k'), hold on,
        xlabel('Rew Diff'), ylabel('p chose higher effort'), title(sprintf('y50 %.2f,x50 %.2f',out.yValueAtSlopePeak,out.fitParams(1)))
        axis square

        if DO_SAVEFIGURES
            currentDateNumber = now;
            currentDateNumber = datestr(currentDateNumber, 'yyyymmdd');
            figurefilename = sprintf('fig_ec_AcrossSessions_type1_%s_%s.pdf',iMonkey,currentDateNumber);
            saveas(gcf, [FIGURE_Folder filesep figurefilename], 'pdf');
        end
        close gcf

    end



end

sprintf('end of script.\n'),













