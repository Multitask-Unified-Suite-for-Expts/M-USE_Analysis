function parsedData = ParseSyncboxData(syncBoxSerialRecvData)

commandStrings = {'HLP', 'ECH', 'IDQ', 'INI', 'LOG', 'LIN', 'LVB', 'TPW', 'TPP', 'TIM', 'RWD', ...
    'NSU', 'NHD', 'NPD', 'NEU', 'NDW', 'CSL', 'CSR', 'CSI', 'CAO', 'CAF', 'CTR', 'CTL', 'FIL', 'FST', 'FSW', 'XAL', 'XJL'};

% parsedData = struct('Commands', {}, 'JoyAndPhoto', {}, 'Sync', {}, 'Reward', {}, 'EventCode', {}, 'SystemStatus', {}, 'Unknown', {});
structData.JoyAndPhotoDetails = struct('Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, ...
    'JoyX', {}, 'JoyY', {}, 'JoyC', {}, 'PhotoL', {}, 'PhotoR', {}, 'NA_PhotoLStatus', {}, 'NA_PhotoRStatus', {});
structData.SyncDetails = struct('Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {});
structData.RewardDetails = struct('Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, 'RewardDuration', {});
structData.Reward2Details = struct('Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, 'RewardDuration', {});
structData.EventCodeDetails = struct('Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, 'EventCode', {});
structData.Commands = struct ('Frame', {}, 'FrameStart', {}, 'Message', {});
structData.SystemStatus = struct ('Frame', {}, 'FrameStart', {}, 'Message', {});
structData.Unknown = struct ('Frame', {}, 'FrameStart', {}, 'Message', {});

parsedData = syncBoxSerialRecvData;

toMatch = [];
reverseStr = '';
for i = 1:height(syncBoxSerialRecvData.Raw)
    
    %print percentage of file parse
    percentDone = 100 * i / height(syncBoxSerialRecvData.Raw);
    msg = sprintf('......Percent done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    %if the previous line was incomplete, this line should be appended to
    %it
    if ~isempty(toMatch)
        toMatch.Message = {[toMatch.Message{:} syncBoxSerialRecvData.Raw.Message{i}]};
        thisLine = toMatch;
    else
        thisLine = syncBoxSerialRecvData.Raw(i,:);
    end
    
    %compare the first three characters of the string to neurarduino 
    %command messages
    message = thisLine.Message{:};
    if length(message) > 2
        if sum(cell2mat(strfind(commandStrings,message(1:3)))) > 0
            structData.Commands = ParseGenericDetails(structData.Commands, thisLine);
        else
            [structData, toMatch] = ParseNonCommandMessages(structData, thisLine);
        end
    elseif ~isempty(message)
        [structData, toMatch] = ParseNonCommandMessages(structData, thisLine);
    end
end

structFields = fields(structData);

for i = 1:length(structFields)
    if ~isempty(structData.(structFields{i}))
        parsedData.(structFields{i}) = struct2table(structData.(structFields{i}));
    else
        parsedData.(structFields{i}) = NaN;
    end
end


% structData.JoyAndPhotoDetails = struct('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, ...
%     'JoyX', {}, 'JoyY', {}, 'JoyC', {}, 'PhotoL', {}, 'PhotoR', {}, 'PhotoLVerdict', {}, 'PhotoRVerdict', {});
% structData.SyncDetails = struct('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {});
% structData.RewardDetails = struct('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, 'RewardDuration', {});
% structData.EventCodeDetails = struct('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'ArduinoTimestamp', {}, 'EventCode', {});
% structData.Commands = struct ('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'Message', {});
% structData.SystemStatus = struct ('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'Message', {});
% structData.Unknown = struct ('Subject', {}, 'Frame', {}, 'FrameStart', {}, 'Message', {});


function [structData, toMatch] = ParseNonCommandMessages(structData, thisLine)
message = thisLine.Message{:};
toMatch = [];

%compare to the possible message types
switch message(1)
    case 'T' %standard neurarduino joystick + photo data line (every 3.3 ms)
        if length(message) == 23 && (strcmp(message(end),'B') || strcmp(message(end),'W'))
            structData.JoyAndPhotoDetails = ParseJoyAndPhotoDetails(structData.JoyAndPhotoDetails, thisLine);
        elseif length(message) < 23
            toMatch = thisLine;
        else
            structData.Unknown = ParseGenericDetails(structData.Unknown, thisLine);
        end
    case 'S' %either the first line of a system status report, or a sync data line
        if length(message) > 2
            if strcmp(message(1:3), 'Sys') && length(message) == 37
                structData.SystemStatus = ParseGenericDetails(structData.SystemStatus, thisLine);
            elseif strcmp(message(1:3), 'Syn') && length(message) == 15
                structData.Sync = ParseSyncDetails(structData.Sync, thisLine);
            elseif length(message) < 15
                toMatch = thisLine;
            else
                structData.Unknown = ParseGenericDetails(structData.Unknown, thisLine);
            end
        else
            toMatch = thisLine;
        end
    case 'R' %reward
        if length(message) == 21 
            structData.RewardDetails = ParseRewardDetails(structData.RewardDetails, thisLine);
        elseif length(message) == 23
            structData.Reward2Details = ParseReward2Details(structData.Reward2Details, thisLine);
        elseif length(message) < 21
            toMatch = thisLine;
        else
            structData.Unknown = ParseGenericDetails(structData.Unknown, thisLine);
        end
    case 'C' %event code ("lynx")
        if length(message) == 19
            structData.EventCodeDetails = ParseEventCodeDetails(structData.EventCodeDetails, thisLine);
        elseif length(message) < 19
            toMatch = thisLine;
        else
            structData.Unknown = ParseGenericDetails(structData.Unknown, thisLine);
        end 
    case ' ' %middle of system status report
        structData.SystemStatus = ParseGenericDetails(structData.SystemStatus, thisLine);
    case 'E' %end of system status report
        structData.SystemStatus = ParseGenericDetails(structData.SystemStatus, thisLine);
    otherwise
        structData.Unknown = ParseGenericDetails(structData.Unknown, thisLine);
end

function structData = ParseGenericDetails(structData, newData)
% structData(1).Subject = [structData.Subject; newData.Subject];
structData(1).Frame = [structData.Frame; newData.FrameRecv];
structData(1).FrameStart = [structData.FrameStart; newData.FrameStart] * 1000;
structData(1).Message = [structData.Message; newData.Message];

function joyAndPhotoDetails = ParseJoyAndPhotoDetails(joyAndPhotoDetails, newData)
% joyAndPhotoDetails(1).Subject = [joyAndPhotoDetails.Subject; newData.Subject];
joyAndPhotoDetails(1).Frame = [joyAndPhotoDetails.Frame; newData.FrameRecv];
joyAndPhotoDetails(1).FrameStart = [joyAndPhotoDetails.FrameStart; newData.FrameStart * 1000];
joyAndPhotoDetails(1).ArduinoTimestamp = [joyAndPhotoDetails.ArduinoTimestamp; hex2dec(newData.Message{1}(2:9)) / 10];
joyAndPhotoDetails(1).JoyX = [joyAndPhotoDetails.JoyX; hex2dec(newData.Message{1}(11:12))];
joyAndPhotoDetails(1).JoyY = [joyAndPhotoDetails.JoyY; hex2dec(newData.Message{1}(13:14))];
joyAndPhotoDetails(1).JoyC = [joyAndPhotoDetails.JoyC; hex2dec(newData.Message{1}(15:16))];
joyAndPhotoDetails(1).PhotoL = [joyAndPhotoDetails.PhotoL; hex2dec(newData.Message{1}(18:19))];
joyAndPhotoDetails(1).PhotoR = [joyAndPhotoDetails.PhotoR; hex2dec(newData.Message{1}(20:21))];
if strcmp(newData.Message{1}(22), 'B')
    joyAndPhotoDetails(1).NA_PhotoLStatus = [joyAndPhotoDetails.NA_PhotoLStatus; 0];
elseif strcmp(newData.Message{1}(22), 'W')
    joyAndPhotoDetails(1).NA_PhotoLStatus = [joyAndPhotoDetails.NA_PhotoLStatus; 1];
end
if strcmp(newData.Message{1}(23), 'B')
    joyAndPhotoDetails(1).NA_PhotoRStatus = [joyAndPhotoDetails.NA_PhotoRStatus; 0];
elseif strcmp(newData.Message{1}(23), 'W')
    joyAndPhotoDetails(1).NA_PhotoRStatus = [joyAndPhotoDetails.NA_PhotoRStatus; 1];
end

function syncDetails = ParseSyncDetails(syncDetails, newData)
% syncDetails(1).Subject = [syncDetails.Subject; newData.Subject];
syncDetails(1).Frame = [syncDetails.Frame; newData.FrameRecv];
syncDetails(1).FrameStart = [syncDetails.FrameStart; newData.FrameStart * 1000];
syncDetails(1).ArduinoTimestamp = [syncDetails.ArduinoTimestamp; hex2dec(newData.Message{1}(8:15)) / 10];

function rewardDetails = ParseRewardDetails(rewardDetails, newData)
% rewardDetails(1).Subject = [rewardDetails.Subject; newData.Subject];
rewardDetails(1).Frame = [rewardDetails.Frame; newData.FrameRecv];
rewardDetails(1).FrameStart = [rewardDetails.FrameStart; newData.FrameStart * 1000];
rewardDetails(1).ArduinoTimestamp = [rewardDetails.ArduinoTimestamp; hex2dec(newData.Message{1}(9:16)) / 10];
rewardDetails(1).RewardDuration = [rewardDetails.RewardDuration; hex2dec(newData.Message{1}(18:21))];

function rewardDetails = ParseReward2Details(rewardDetails, newData)
% rewardDetails(1).Subject = [rewardDetails.Subject; newData.Subject];
rewardDetails(1).Frame = [rewardDetails.Frame; newData.FrameRecv];
rewardDetails(1).FrameStart = [rewardDetails.FrameStart; newData.FrameStart * 1000];
rewardDetails(1).ArduinoTimestamp = [rewardDetails.ArduinoTimestamp; hex2dec(newData.Message{1}(11:18)) / 10];
rewardDetails(1).RewardDuration = [rewardDetails.RewardDuration; hex2dec(newData.Message{1}(20:23))];

function eventCodeDetails = ParseEventCodeDetails(eventCodeDetails, newData)
% eventCodeDetails(1).Subject = [eventCodeDetails.Subject; newData.Subject];
eventCodeDetails(1).Frame = [eventCodeDetails.Frame; newData.FrameRecv];
eventCodeDetails(1).FrameStart = [eventCodeDetails.FrameStart; newData.FrameStart * 1000];
eventCodeDetails(1).ArduinoTimestamp = [eventCodeDetails.ArduinoTimestamp; hex2dec(newData.Message{1}(7:14)) / 10];
eventCodeDetails(1).EventCode = [eventCodeDetails.EventCode; hex2dec(newData.Message{1}(16:19))];




