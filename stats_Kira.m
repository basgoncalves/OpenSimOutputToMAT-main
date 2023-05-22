
function stats_Kira()


add_repos_to_path
[S,Results] = get_data_struct;

Vars = {'IK','ID','SO','SO_Activation','JRL'};

for iSubj = 1:length(S.subjects)
    for iSess = 1:length(S.sessions)
        
        Paths = get_data_paths(iSubj,iSess);
        if ~isfolder(Paths.session)
            disp([Paths.session ' does not exist'])
            continue
        end
        
        load(Paths.matResults)
         
        for iVar = 1:length(Vars)
            curr_var = Vars{iVar};
            data_struct = data.(curr_var).FINAL_PERSONALISEDTORSIONS_scaled_final;
            results_struct = Results.(curr_var).(['session' num2str(iSess)]);
            [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct);
            
            Results.(curr_var).(['session' num2str(iSess)]) = results_struct;

        end       
    end
end

Results
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
leg = {'right','left'};
for i = 1:3
    for l = 1:numel(leg)
        Results.IK.(['session' num2str(i)]).(leg{l}) = struct;
        Results.ID.(['session' num2str(i)]).(leg{l}) = struct;
        Results.SO.(['session' num2str(i)]).(leg{l}) = struct;
        Results.SO_Activation.(['session' num2str(i)]).(leg{l}) = struct;
        Results.JRL.(['session' num2str(i)]).(leg{l}) = struct;
    end
end


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
function [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct)
% data_struct = a struct similar to the output of load_sto_file

trialNames = fields(data_struct);
if contains(trialNames{1},'mass')
    trialNames(1) = [];
end
variables  = fields(data_struct.(trialNames{1}));
variables = variables(~contains(variables,'calc')); % remove variables that contain calc


raw   = struct; raw.left = struct; raw.right = struct;
tnorm = struct; tnorm.left = struct; tnorm.right = struct;


% crete the full struct;
for iVar = variables'
    raw.right.(iVar{1}) = [];
    raw.left.(iVar{1}) = [];
    tnorm.right.(iVar{1}) = [];
    tnorm.left.(iVar{1}) = [];
end


% separate right and left trials into different structs
for iTrial = trialNames'
    for iVar = variables'
        try
            currentTrialData = data_struct.(iTrial{1}).(iVar{1});  % if variable doesn't exist in one trial, make it NaN
        catch
            currentTrialData = nan(101,1);
        end
        rows = length(currentTrialData);
        t1 = data_struct.(iTrial{1}).time(1);
        t2 = data_struct.(iTrial{1}).time(2);
        fs = 1/(t2-t1);
        if contains(iTrial{1},'_right')
            raw.right.(iVar{1})(1:rows,end+1) = currentTrialData;
            tnorm.right.(iVar{1})(1:101,end+1) = TimeNorm(currentTrialData,fs);
        else
            raw.left.(iVar{1})(1:rows,end+1) = currentTrialData;
            tnorm.left.(iVar{1})(1:101,end+1) = TimeNorm(currentTrialData,fs);
        end
    end
end

% make zeros = NaNfor iTrial = trialNames'
for iVar = variables'
    raw.right.(iVar{1})(raw.right.(iVar{1})==0) = NaN;
    raw.left.(iVar{1})(raw.left.(iVar{1})==0) = NaN;
end


% mean tnorm data 
for iVar = variables'
    tnorm.right.(iVar{1}) = mean(tnorm.right.(iVar{1}),2);
    tnorm.left.(iVar{1})  = mean(tnorm.left.(iVar{1}),2);
end


% add to results struct
for iVar = variables'
    if ~isfield(results_struct.right,iVar{1})
        results_struct.right.(iVar{1}) = [];
        results_struct.left.(iVar{1}) = [];
    end
    results_struct.right.(iVar{1})(:,end+1) = tnorm.right.(iVar{1});
    results_struct.left.(iVar{1})(:,end+1) = tnorm.left.(iVar{1});
end


