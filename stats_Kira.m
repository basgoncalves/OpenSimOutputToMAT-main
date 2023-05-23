
function stats_Kira()


add_repos_to_path
[Results] = get_data_struct; 
S = get_subjects;

gather_data_in_struct = true;
run_stats = true;

%% organise data in struct 
if gather_data_in_struct
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
    cd(get_main_dir)
    save('results.mat', 'Results')
end

%% stats
if run_stats
    load([get_main_dir fp 'results.mat'])
    
    Results = calculate_peaks(Results);
    load_participant_characteristics()

    x = Results.JRL.session1.left.peak_HCF';
    y = Results.JRL.session1.right.peak_HCF';
    [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,1,0.05);

    
  

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
function data_path = get_main_dir()
activeFile = [mfilename('fullpath') '.m'];
data_path = [fileparts(fileparts(activeFile)) '\Kira_MSc_data'];

% --------------------------------------------------------------------------------------------------------------- %
function S = get_subjects()
S = struct;
S.subjects = {getfolders(get_main_dir).name};
S.sessions = {'\pre', '\post', ''};

% --------------------------------------------------------------------------------------------------------------- %
function load_participant_characteristics()
cd(get_main_dir)
load('participants_characteristics.mat')

% --------------------------------------------------------------------------------------------------------------- %
function [Results] = get_data_struct()

Results = struct;
leg = {'right','left'};
for i = 1:3
    for l = 1:numel(leg)
        session = ['session' num2str(i)];
        Results.IK.(session).(leg{l}) = struct;
        Results.ID.(session).(leg{l}) = struct;
        Results.SO.(session).(leg{l}) = struct;
        Results.SO_Activation.(session).(leg{l}) = struct;
        Results.JRL.(session).(leg{l}) = struct;
        Results.Mass.(session) = [];
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function Paths = get_data_paths(iSubj,iSess)

S = get_subjects;
session_path = [get_main_dir fp S.subjects{iSubj} fp S.sessions{iSess} fp 'output automization\FINAL_PERSONALISEDTORSIONS_scaled_final'];
trialNames      = {getfolders(session_path).name};

Paths.main      = get_main_dir;
Paths.session   = session_path;
Paths.trials    = cellfun(@(c)[session_path fp c],trialNames,'uni',false);
Paths.trialNames = trialNames;
Paths.matResults = [fileparts(session_path) fp 'dataStruct_ErrorScores_no_trials_removed.mat'];

% --------------------------------------------------------------------------------------------------------------- %
function [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct)
% data_struct = a struct similar to the output of load_sto_file


% add mass to results_struct AND remove "

trialNames = fields(data_struct);
if contains(trialNames{1},'mass')
    trialNames(1) = [];
end

% remove variables that contain calc
variables  = fields(data_struct.(trialNames{1}));
variables = variables(~contains(variables,'calc')); 

% create structs
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
    
    % right leg (try in case the "iVar" is empty)
    try
        results_struct.right.(iVar{1})(:,end+1) = tnorm.right.(iVar{1});
    end
    % leg leg
    try
        results_struct.left.(iVar{1})(:,end+1) = tnorm.left.(iVar{1});
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function Results = calculate_peaks(Results)

Legs = {'right','left'};

for i = 1:3
    for l = 1:length(Legs)
        session = ['session' num2str(i)];
        leg = Legs{l};

        % find peak HCF
        Results.JRL.(session).(leg).HCF = sum3D(Results.JRL.(session).(leg), ['hip_' leg(1) '_on_pelvis_in_pelvis']);
        Results.JRL.(session).(leg).peak_HCF = calc_max(Results.JRL.(session).(leg).HCF);
        
        % find peak KCF
        Results.JRL.(session).(leg).KCF = sum3D(Results.JRL.(session).(leg), ['knee_' leg(1) '_on_tibia_' leg(1) '_in_tibia_' leg(1)]);
        Results.JRL.(session).(leg).peak_KCF = calc_max(Results.JRL.(session).(leg).KCF);

        Muscles = fields(Results.SO.(session).(leg));
        for m = 1:length(Muscles)
            muscle = Muscles{m};
            Results.SO.(session).(leg).(['peak_' muscle]) = calc_max(Results.SO.(session).(leg).(muscle));
        end

    end
end

% --------------------------------------------------------------------------------------------------------------- %
function out = sum3D(Results_substruct, variable)
x = Results_substruct.([variable '_fx']);
y = Results_substruct.([variable '_fy']);
z = Results_substruct.([variable '_fz']);

out = sqrt(x.^2 + y.^2 + z.^2);

%------------------ calc max moving average 
function out = calc_max(time_curve)

out = max(movmean(time_curve,5));


