
function stats_Kira()


add_repos_to_path
[S,Results] = get_data_struct;


for iSubj = 1:length(S.subjects)
    for iSess = 1:length(S.sessions)
        
        Paths = get_data_paths(iSubj,iSess);
        if ~isfolder(Paths.session)
            disp([Paths.session ' does not exist'])
            continue
        end
        
        load(Paths.matResults)
        [right,left] = mean_left_right_legs(data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final);

        
%         .T_Dynamic08_1_right.ankle_angle_r);
        
    end
end


%% -------------------------------------------------------------------------------------------------------------- %
% ---------------------------------------------------- FUCNTIONS ------------------------------------------------ %
% --------------------------------------------------------------------------------------------------------------- %
function add_repos_to_path
activeFile = [mfilename('fullpath') '.m'];

OpenSimOutputToMAT_Dir = fileparts(activeFile);
addpath(genpath(OpenSimOutputToMAT_Dir))
disp([OpenSimOutputToMAT_Dir ' added to path'])


msk_modellingDir  = [fileparts(OpenSimOutputToMAT_Dir) '\MSKmodelling'];
addpath(genpath(msk_modellingDir))
disp([msk_modellingDir ' added to path'])

% --------------------------------------------------------------------------------------------------------------- %
function [S,Results] = get_data_struct()

activeFile = [mfilename('fullpath') '.m'];

S = struct;
S.dataDir =  [fileparts(fileparts(activeFile)) '\Kira_MSc_data'];
S.subjects = {getfolders(S.dataDir).name};
S.sessions = {'\pre', '\post', ''};

Results = struct;

% --------------------------------------------------------------------------------------------------------------- %
function Paths = get_data_paths(iSubj,iSess)

[S,Results] = get_data_struct();

session_path = [S.dataDir fp S.subjects{iSubj} fp S.sessions{iSess} fp 'output automization\FINAL_PERSONALISEDTORSIONS_scaled_final'];

Paths.session   = session_path;
trialNames      = {getfolders(session_path).name};
Paths.trials    = cellfun(@(c)[session_path fp c],trialNames,'uni',false);
Paths.trialNames = trialNames;
Paths.matResults = [fileparts(session_path) fp 'dataStruct_ErrorScores_no_trials_removed_issue_fixed.mat'];



% --------------------------------------------------------------------------------------------------------------- %
% function save_log()

% --------------------------------------------------------------------------------------------------------------- %
function [right,left] = mean_left_right_legs(Struct_with_trials)

trialNames = fields(Struct_with_trials);
variables  = fields(Struct_with_trials.(trialNames{1}));

left = struct; right = struct;

% crete the full struct;
for iTrial = trialNames'
    for iVar = variables'
        right.(iVar{1}) = [];
        left.(iVar{1}) = [];
    end
end


% separate right and left trials into different structs
for iTrial = trialNames'
    for iVar = variables'
        currentTrialData = Struct_with_trials.(iTrial{1}).(iVar{1});
        rows = length(currentTrialData);

        if contains(iTrial{1},'_right')
            right.(iVar{1})(1:rows,end+1) = currentTrialData;
        else
            left.(iVar{1})(1:rows,end+1) = currentTrialData;
        end
    end
    right.(iVar{1})(right.(iVar{1})==0) = NaN;     left.(iVar{1})(left.(iVar{1})==0) = NaN;
end


% make zeros = NaN
for iTrial = trialNames'
    for iVar = variables'
          right.(iVar{1})(right.(iVar{1})==0) = NaN; 
    left.(iVar{1})(left.(iVar{1})==0) = NaN;
    end
end



