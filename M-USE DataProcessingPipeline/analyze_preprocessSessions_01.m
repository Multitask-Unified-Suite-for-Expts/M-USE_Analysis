
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
%
% ---
% --- this script prpocesses Data from Session folders and saves trialdata,
% --- blockdata and framedata from each task in a matfile
% ---

% SessionSets = {'Frey-1' };
  %SessionSets = {'Frey-2' }; % ok, finished - last run Jan 29
% SessionSets = {'Frey-3' }; % new since jan 2024 - last run Jan 29

% SessionSets = {'Wotan-1' };
 SessionSets = {'Wotan-2' }; %ongoing since august 2023 - last run Jan 29


% SessionSets = {'Jotun-1' }; Training
% SessionSets = {'Jotun-2' }; % VS_WWW_FL_VS_WWW ...- last run Jan 29 with some initial variable FrameEventCodes

% SessionSets = {'Kyrre-1' }; Training
% SessionSets = {'Kyrre-2' }; %VS_WWW_FL_VS_WWW


% SessionSets = {'Bard-1' };  % the WWW task with WM and effort control it-  it ran feb 1st
% SessionSets = {'Sindri-1' };  % the WWW task with WM and effort control -  it ran feb 1st




doReadFrameData = 0;
% --- specify the session folder and the sessions to analyze
for iO = 1:length(SessionSets)
    FOLDER_SESSION = [];
    if strcmp(SessionSets{iO},'Frey-1')
        EXP_ID = 'VS_EC_FL_MRT';
        Subject = 'Frey';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Frey/VS_EC_FL_MZG';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_VS_EC_FL_MZG_Frey']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
        % Frey signficnatly improived on the guided mazw on August31st and Sept st 2023 (after three weeks)!
    elseif strcmp(SessionSets{iO},'Frey-2')
        EXP_ID = 'WM_EC_FL_MZG';
        Subject = 'Frey';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Frey/WM_EC_FL_MZG';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_WM_EC_FL_MZG_Frey']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
elseif strcmp(SessionSets{iO},'Frey-3')
        EXP_ID = 'MRT_FL_VS_WM';
        Subject = 'Frey';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Frey/MRT_FL_VS_WM';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_MRT_FL_VS_WM_Frey']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end

    elseif strcmp(SessionSets{iO},'Wotan-1')
        EXP_ID = 'VS_EC_FL_MRT';
        Subject = 'Wotan';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Wotan/VS_EC_FL_MZG';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_VS_EC_FL_MRT_Wotan']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
        % FOLDER_SESSION{iSession}  = 'Session_08_14_23__09_17_16_Wotan'; % has errorneous blockdata
        % FOLDER_SESSION{iSession}  = 'Session_08_15_23__09_36_55_Wotan'; % first with TrialID and with good ID.ED, old FL block sequence
        % FOLDER_SESSION{iSession}  = 'Session_08_17_23__09_53_45_Wotan'; % new FL blocksequence (60% EDsame/IDsame) blocks
    elseif strcmp(SessionSets{iO},'Wotan-2')
        Subject = 'Wotan';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Wotan/WM_EC_FL_MZG';
        EXP_ID = 'WM_EC_FL_MZG';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_WM_EC_FL_MZG_Wotan']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
        % FOLDER_DATA = /Volumes/'Womelsdorf Lab'/DATA_kiosk/Wotan/WM_EC_FL_MZG/
        % FOLDER_SESSION{iSession}  = 'Session_08_14_23__09_17_16_Wotan'; % has errorneous blockdata
        % FOLDER_SESSION{iSession}  = 'Session_08_15_23__09_36_55_Wotan'; % first with TrialID and with good ID.ED, old FL block sequence
        % FOLDER_SESSION{iSession}  = 'Session_08_17_23__09_53_45_Wotan'; % new FL blocksequence (60% EDsame/IDsame) blocks
        % FOLDER_SESSION{iSession}  = 'Session_08_25_23__08_35_30_Wotan'; % new FL blocksequence (60% EDsame/IDsame) blocks
        % --- collect mat files
    elseif strcmp(SessionSets{iO},'Jotun-1')
        EXP_ID = 'VS_FL';
        Subject = 'Jotun';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Jotun/VS_FL';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_VS_FL_Jotun']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
    elseif strcmp(SessionSets{iO},'Jotun-2')
    
        EXP_ID = 'VS_WWW_FL_VS_WWW';
        Subject = 'Jotun';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Jotun/VS_WWW_FL_VS_WWW';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_VS_WWW_FL_VS_WWW_Jotun']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
    
    elseif strcmp(SessionSets{iO},'Kyrre-1')
        EXP_ID = 'VS_FL';
        Subject = 'Kyrre';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Kyrre/VS_FL';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_VS_FL_Kyrre']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
    elseif strcmp(SessionSets{iO},'Kyrre-2')
        EXP_ID = 'VS_WWW_FL_VS_WWW';
        Subject = 'Kyrre';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Kyrre/VS_WWW_FL_VS_WWW';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_VS_WWW_FL_VS_WWW_Kyrre']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end

    elseif strcmp(SessionSets{iO},'Bard-1')
        Subject = 'Bard';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Bard/MUSE/New_WWW_WM';
        EXP_ID = 'WWW_WM_WWW';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_WWW_WM_WWW_Bard']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
    
    elseif strcmp(SessionSets{iO},'Sindri-1')
   Subject = 'Sindri';
        FOLDER_DATA = '/Volumes/Womelsdorf Lab/DATA_kiosk/Sindri/MUSE/New_WWW_WM';
        EXP_ID = 'WWW_WM_WWW';
        RESULT_FOLDER = [pwd filesep 'MUSEMAT01_WWW_WM_WWW_Sindri']; if ~exist(RESULT_FOLDER), mkdir(RESULT_FOLDER), end
    
    end

    % --- --- --- --- --- --- ---
    % --- collect all files in the Datafolder if they were not specified
    % --- --- --- --- --- --- ---
    if isempty(FOLDER_SESSION)
        dirInfo = dir(FOLDER_DATA);
        FOLDER_SESSION = {dirInfo([dirInfo.isdir]).name};
        if isempty(dirInfo), sprintf('could not open %s (FOLDER_DATA)\n',FOLDER_DATA), continue, end
        FOLDER_SESSION = FOLDER_SESSION(~ismember(FOLDER_SESSION, {'.', '..'}));
        % --- for processing the latest data first...
        FOLDER_SESSION  = fliplr(FOLDER_SESSION);
    end

    % --- --- --- --- --- --- ---
    % --- determine which sessions to analyze
    % --- --- --- --- --- --- ---
    for iSession = 56:85 % EDITED FOR SESSION FROM 2/8 TO 3/21
        iSessionName = FOLDER_SESSION{iSession};
        iSessionDataFolder = [ FOLDER_DATA  filesep  iSessionName ];
        iResultFile = ['DAT01_' Subject '_' EXP_ID '_' iSessionName];
        iResultFileWithPath = [RESULT_FOLDER filesep iResultFile];

        % --- --- --- --- --- --- --- --- --- --- ---
        % --- pre-process M-USE session folder with all tasks
        % --- --- --- --- --- --- --- --- --- --- ---
        iDataFolder = [ iSessionDataFolder ];
        if ~exist([iSessionDataFolder filesep 'ProcessedData/'])
            Import_MUSE_Session_RawData('dataFolder', [iDataFolder],'gazeArgs','cancel','serialDataArgs','cancel');
        end

        % --- --- --- --- --- --- ---
        % --- find the tasks available in this session folder
        % --- --- --- --- --- --- ---
        taskNames = {};
        %
        dirInfo = dir(iSessionDataFolder)
        folderNames = {dirInfo([dirInfo.isdir]).name};
        if isempty(dirInfo), sprintf('could not open %s\n',iDataFolderSession), continue, end
        folderNames = folderNames(~ismember(folderNames, {'.', '..'})); % Removes '.' (current directory) and '..' (parent directory) from the list
        pattern = '^Task\d{4}_\w+';% Regular expression pattern to match TaskXXXX_Name
        taskNames = regexp(folderNames, pattern, 'match');% Use cellfun and regexp to extract task names
        taskNames = vertcat(taskNames{:});  % Convert cell array of cell arrays to a single cell array

        disp('...collecting data of these tasks:');
        disp(taskNames);

        if isempty(taskNames), 
            disp('empty taskNames, returning'), return
        end
        
        dat = [];
        dat.sessionFolder  = iSessionDataFolder;
        dat.sessionName    = iSessionName;
        dat.sessionTasks   = taskNames;
        dat.subject         = Subject;
        dat.expID           = EXP_ID;

        dat.taskLabel =[];
        dat.trialData = [];
        dat.blockData = [];
        dat.frameData = [];

        dat.cfg_BlockDef = [];
        dat.cfg_StimDef = [];
        dat.cfg_TrialDef = [];
        dat.cfg_MazeDef = [];

        % --- --- --- --- --- --- --- --- --- --- ---
        % --- read  trial, block and framdata
        % --- --- --- --- --- --- --- --- --- --- ---
        for iT = 1:length(taskNames)

            dat.taskLabel{iT} = taskNames{iT};

            % --- --- --- --- --- --- --- --- --- --- ---
            % --- read  trial, block and framdata
            % --- --- --- --- --- --- --- --- --- --- ---
            preFixFolder = [iSessionDataFolder filesep 'ProcessedData' filesep taskNames{iT} filesep];

            GotData = 0;
            if exist([preFixFolder  'TrialData.mat'])
                in = load([preFixFolder  'TrialData.mat']);
                dat.trialData{iT} = table2struct(in.trialData);
                GotData=1;
            end

            if exist([preFixFolder  'BlockData.mat'])
                in = load([preFixFolder  'BlockData']);
                dat.blockData{iT}  = table2struct(in.blockData);
                GotData=1;
            end

            if doReadFrameData == 1 & exist([preFixFolder  'FrameData'])~=0
                in = load([preFixFolder   'FrameData']);
                dat.frameData{iT}  = table2struct(in.frameData);
            end

            % --- do not add information when neither trial nor block data was found.
            if GotData == 0, continue, end

            % --- --- --- --- --- --- --- --- --- --- ---
            % --- read session configs for each task too.
            % --- --- --- --- --- --- --- --- --- --- ---
            configFolderInfo = dir([iSessionDataFolder filesep  'SessionConfigs' filesep]);
            folderNames = {configFolderInfo([configFolderInfo.isdir]).name};
            folderNames = folderNames(~ismember(folderNames, {'.', '..'}))'; % Removes '.' (current directory) and '..' (parent directory) from the list
            for iConfig=1:length(folderNames), if contains(taskNames{iT}, folderNames{iConfig})==1, break, end, end

            configFiles = dir([iSessionDataFolder filesep  'SessionConfigs' filesep folderNames{iConfig}]);
            for j=1:length(configFiles)
                if (configFiles(j).isdir) | strcmp(configFiles(j).name(j),'.' ),  continue, end
                if ~isempty(findstr(configFiles(j).name,'BlockDef'))
                    dat.cfg_BlockDef{iT} = readtable([configFiles(j).folder filesep configFiles(j).name], 'Delimiter', '\t');
                elseif ~isempty(findstr(configFiles(j).name,'StimDef'))
                    dat.cfg_StimDef{iT} = readtable([configFiles(j).folder filesep configFiles(j).name], 'Delimiter', '\t');
                elseif ~isempty(findstr(configFiles(j).name,'TrialDef'))
                    dat.cfg_TrialDef{iT} = readtable([configFiles(j).folder filesep configFiles(j).name], 'Delimiter', '\t');
                elseif ~isempty(findstr(configFiles(j).name,'MazeDef'))
                    dat.cfg_MazeDef{iT} = readtable([configFiles(j).folder filesep configFiles(j).name], 'Delimiter', '\t');
                end
            end
        end

%dat.taskLabel

        % --- --- --- --- --- --- --- --- --- --- --- --- ---
        % --- save the preprocessed session data
        % --- --- --- --- --- --- --- --- --- --- --- --- ---
        save(iResultFileWithPath,'dat','-V7.3')
        sprintf('preprocessed sesssion %d of %d; saved dat in %s\n',iSession, length(FOLDER_SESSION), iResultFile),

    end

    sprintf('done with %d sessions.\n',length(iSession)),
end

return

