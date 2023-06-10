function stats_Kira()
addpath(genpath('C:\Users\Balu\Nextcloud\Documents\MA\Code\MSKmodelling'))

% add_repos_to_path
[Results] = get_data_struct;
S = get_subjects;

gather_data_in_struct = false;

run_stats =                 true;
plot_JRF_curve_and_peaks =  false;
scatter_peaks_angles =      false; two_peaks = false;
plot_corr =                 false;
multiple_regress =          false;

%% organise data in struct
if gather_data_in_struct
    Analyses = {'IK','ID','SO','SO_Activation','JRL'};
    for iSubj = 1:length(S.subjects)
        for iSess = 1:length(S.sessions)

            Paths = get_data_paths(iSubj,iSess);
            if ~isfolder(Paths.session)
                disp([Paths.session ' does not exist'])
                continue
            end

            load(Paths.matResults)

            for iVar = 1:length(Analyses)
                curr_analysis = Analyses{iVar};
                data_struct = data.(curr_analysis).FINAL_PERSONALISEDTORSIONS_scaled_final;
                results_struct = Results.(curr_analysis).(['session' num2str(iSess)]);

                [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct);

                Results.(curr_analysis).(['session' num2str(iSess)]) = results_struct;

            end
        end
    end
    cd(get_main_dir)
    save('results.mat', 'Results')
end

%% stats
if run_stats
    load([get_main_dir fp 'results.mat'])

    legs = {'left'; 'right'};
sessions = {'session1'; 'session2'; 'session3'};
    joints = {'HCF';'KCF';'ACF'};
        angles = {'NSA'; 'AVA'; 'TT'};
    lgnd = {'P01', 'P02', 'P03', 'P04', 'P05';'P01', 'P02', 'P03', 'P04', 'P05';'TD02', 'TD03', 'TD04','',''};
   

    subj_charac = importdata('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\participants_characteristics.mat');
    add_nan = {'ercspn_r'; 'ercspn_l'; 'intobl_r'; 'intobl_l'; 'extobl_r'; 'extobl_l'};
    for k = 1:size(legs,1)
        for i = 1:size(add_nan,1)
            for c=4:-1:2
                Results.SO.session3.(legs{k}).(add_nan{i})(:,c) =Results.SO.session3.(legs{k}).(add_nan{i})(:,c-1);
            end
            Results.SO.session3.(legs{k}).(add_nan{i})(:,1) = NaN;
        end
    end

    Results = calculate_peaks(Results);

    

    % normalize results to bodyweight
    Results_BW = struct;
    for i=1:size(sessions,1)
        for j=1:size(legs,1)
            JRL_fields = fields(Results.JRL.(sessions{i}).(legs{j}));
            for k=1:size(JRL_fields,1)
                if contains(JRL_fields{k},'_loc')
                    Results_BW.JRL.(sessions{i}).(legs{j}).(JRL_fields{k}) = Results.JRL.(sessions{i}).(legs{j}).(JRL_fields{k});
                else
                    Results_BW.JRL.(sessions{i}).(legs{j}).(JRL_fields{k}) = (Results.JRL.(sessions{i}).(legs{j}).(JRL_fields{k}))./(subj_charac.(sessions{i}).mass*9.81);
                end
            end
            SO_fields = fields(Results.SO.(sessions{i}).(legs{j}));

            for k = 1:size(SO_fields,1)
                if contains(SO_fields{k},'_loc')
                    Results_BW.SO.(sessions{i}).(legs{j}).(SO_fields{k}) = Results.SO.(sessions{i}).(legs{j}).(SO_fields{k});
                else
                    Results_BW.SO.(sessions{i}).(legs{j}).(SO_fields{k}) = (Results.SO.(sessions{i}).(legs{j}).(SO_fields{k}))./(subj_charac.(sessions{i}).mass*9.81);
                end
            end
        end
    end

    % plot JRF by session incl. peaks
    if plot_JRF_curve_and_peaks
        nb_plots = size(sessions,1)*size(joints,1);
        for j=1:size(legs,1)
            [ax, pos,FirstCol,LastRow,LastCol] = tight_subplot(nb_plots,0);
            count = 1;
            for iDOF = 1:size(sessions,1)
                for k=1:size(joints,1)
                                axes(ax(count)); hold on; grid on;
                    plot(Results_BW.JRL.(sessions{iDOF}).(legs{j}).(joints{k})); hold on;
                    plot(Results_BW.JRL.(sessions{iDOF}).(legs{j}).(['peak_' joints{k} '_loc']), Results_BW.JRL.(sessions{iDOF}).(legs{j}).(['peak_' joints{k} '_val']),'*')
                    title([sessions{iDOF} ' ' legs{j} ' ' joints{k}])
                    legend(lgnd{iDOF,:})
                    xlim([0, 100])
                    if any(count == FirstCol)
                        ylabel('Joint contact force [N/BW]')
                    elseif any(count == LastRow)
                        xlabel('Gait Cycle [%]')
                    end
                    count = count +1;
                end
            end
                tight_subplot_ticks (ax,LastRow,0)
        end
    end

    % scatter plots peaks and AVA/NSA/TT
    if scatter_peaks_angles
    bone_geom = struct;
    for i = 1:size(angles,1)
        bone_geom.(angles{i}) = [];
        for k =1: size(sessions,1)
            bone_geom_1session = [subj_charac.(sessions{k}).(legs{1}).(angles{i}), subj_charac.(sessions{k}).(legs{2}).(angles{i})];
            bone_geom.(angles{i}) = [bone_geom.(angles{i}), bone_geom_1session];
        end
    end

    nb_plots = size(angles,1)*size(joints,1);
    [ax, pos,FirstCol,LastRow,LastCol] = tight_subplot(nb_plots,0,[],[0.1 0.05],[0.1 0.05]);
    count = 1;
    for j=1:size(joints,1)
        for k=1:size(angles,1)
            axes(ax(count)); hold on; grid on;
                            plot(Results_BW.JRL.(sessions{iDOF}).(legs{j}).(joints{k})); hold on;
                            scatter(bone_geom.(angles{k}) (1:10), [Results_BW.JRL.session1.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session1.(legs{2}).(['peak_' joints{j} '_val'])(1,:)],'r')
                            scatter(bone_geom.(angles{k}) (11:20), [Results_BW.JRL.session2.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session2.(legs{2}).(['peak_' joints{j} '_val'])(1,:)],'b')
                            scatter(bone_geom.(angles{k}) (21:end), [Results_BW.JRL.session3.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session3.(legs{2}).(['peak_' joints{j} '_val'])(1,:)],'k')
                            
                            if two_peaks
                            scatter(bone_geom.(angles{k}) (1:10), [Results_BW.JRL.session1.(legs{1}).(['peak_' joints{j} '_val'])(2,:), Results_BW.JRL.session1.(legs{2}).(['peak_' joints{j} '_val'])(2,:)],'r')
                            scatter(bone_geom.(angles{k}) (11:20), [Results_BW.JRL.session2.(legs{1}).(['peak_' joints{j} '_val'])(2,:), Results_BW.JRL.session2.(legs{2}).(['peak_' joints{j} '_val'])(2,:)],'b')
                            scatter(bone_geom.(angles{k}) (21:end), [Results_BW.JRL.session3.(legs{1}).(['peak_' joints{j} '_val'])(2,:), Results_BW.JRL.session3.(legs{2}).(['peak_' joints{j} '_val'])(2,:)],'k')
                            end

            %                 title([sessions{iDOF} ' ' legs{j} ' ' joints{k}])
            if two_peaks
            legend('pre', '', 'post','', 'TD', '') % two peaks
            else
            legend('pre', 'post', 'TD') % one peak
            end

            if any(count == FirstCol)
                ylabel(['Peak ' joints{j} ' [N/BW]'])
            end
            if any(count == LastRow)
                xlabel([angles{k} ' [°]'])
            end
            count = count +1;
        end
    end
    tight_subplot_ticks (ax,LastRow,1)
    end

    % correlation plots
    if plot_corr
    bone_geom = struct;
    for i = 1:size(angles,1)
        bone_geom.(angles{i}) = [];
        for k =1: size(sessions,1)
            bone_geom_1session = [subj_charac.(sessions{k}).(legs{1}).(angles{i}), subj_charac.(sessions{k}).(legs{2}).(angles{i})];
            bone_geom.(angles{i}) = [bone_geom.(angles{i}), bone_geom_1session];
        end
    end

    nb_plots = size(angles,1)*size(joints,1);
    [ax, pos,FirstCol,LastRow,LastCol] = tight_subplotBG(nb_plots,0,[],[0.1 0.05],[0.1 0.05]);
    count = 1;
    for j=1:size(joints,1)
        for k=1:size(angles,1)
            axes(ax(count)); hold on; grid on;
               scatter(bone_geom.(angles{k}) (21:end), [Results_BW.JRL.session3.(legs{1}).(['peak_' joints{j} '_val'])(2,:), Results_BW.JRL.session3.(legs{2}).(['peak_' joints{j} '_val'])(2,:)],'k')
            y = [Results_BW.JRL.session1.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session1.(legs{2}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session2.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session2.(legs{2}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session3.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session3.(legs{2}).(['peak_' joints{j} '_val'])(1,:)];

            plotCorr(bone_geom.(angles{k})',y',1,0.05)

            %                 title([sessions{iDOF} ' ' legs{j} ' ' joints{k}])
            % legend('pre', '', 'post','', 'TD', '') % two peaks
            % legend('pre', 'post', 'TD') % one peak
            if any(count == FirstCol)
                ylabel(['Peak ' joints{j} ' [N/BW]'])
            end
            if any(count == LastRow)
                xlabel([angles{k} ' [°]'])
            end
            count = count +1;
        end
    end
    tight_subplot_ticks (ax,LastRow,1)
    end

    %     x = Results.JRL.session1.left.peak_HCF';
    %     y = Results.JRL.session1.right.peak_HCF';
    %     [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,1,0.05);


    % plot kinematics
    plotKinematics(Results, 'Absolute')

% multiple regression analysis
if multiple_regress
    % b: coefficient estimates, bint: 95% confidence intervals, r: residuals, rint: intervals to diagonose outliers, stats: R², F, p, error variance
%     X(:,1) = ones(28,1);
    for i = 1:size(angles,1)
        bone_angles = [];
        for k =1: size(sessions,1)
            bone_angles_1session = [subj_charac.(sessions{k}).(legs{1}).(angles{i})'; subj_charac.(sessions{k}).(legs{2}).(angles{i})'];
            bone_angles = [bone_angles; bone_angles_1session];
        end
        X(:,i) = bone_angles;
    end

    for j =1:size(joints,1)
        y{j} = [Results_BW.JRL.session1.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session1.(legs{2}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session2.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session2.(legs{2}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session3.(legs{1}).(['peak_' joints{j} '_val'])(1,:), Results_BW.JRL.session3.(legs{2}).(['peak_' joints{j} '_val'])(1,:)]';
        modelfun = @(b,x)b(1) + b(2)*x(:,1).^b(3) + b(4)*x(:,2).^b(5)+b(6)*x(:,3).^b(7);
        beta0 = [-50 20 -1 150 -1 4 37];
        mdl{j} = fitnlm(X,y{j},modelfun,beta0);
        clear modelfun beta0
    end

% plot data and regression model
% for i = size(y,2)
% scatter3(X(:,1),X(:,2),X(:,3),y{1,i},'filled')
% hold on
% x1fit = min(X(:,1)):100:max(X(:,1));
% x2fit = min(X(:,2)):10:max(X(:,2));
% x3fit = min(X(:,3)):10:max(X(:,3));
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% YFIT = b{1,i}(1) + b{1,i}(2)*X1FIT + b{1,i}(3)*X2FIT + b{1,i}(4)*X1FIT.*X2FIT;
% mesh(X1FIT,X2FIT,YFIT)
% xlabel('Weight')
% ylabel('Horsepower')
% zlabel('MPG')
% view(50,10)
% hold off
% clear X1FIT X2FIT YFIT

end

    
%     Results = calculate_peaks(Results);                                                                             

    load_participant_characteristics()

%     x = Results.JRL.session1.left.peak_HCF';
%     y = Results.JRL.session1.right.peak_HCF';
%     [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,1,0.05);  
% 

end

%% -------------------------------------------------------------------------------------------------------------- %
% ---------------------------------------------------- FUCNTIONS ------------------------------------------------ %
% --------------------------------------------------------------------------------------------------------------- %

function out  = fp
out  = filesep;

% --------------------------------------------------------------------------------------------------------------- %
function add_repos_to_path
activeFile = [mfilename('fullpath') '.m'];

OpenSimOutputToMAT_Dir = fileparts(activeFile);
addpath(genpath(OpenSimOutputToMAT_Dir))
disp([OpenSimOutputToMAT_Dir ' added to path'])


msk_modellingDir  = [fileparts(OpenSimOutputToMAT_Dir) '\MSKmodelling'];
rmpath(genpath(msk_modellingDir))

rmpath(genpath(msk_modellingDir))

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


%         Results.IK.(['session' num2str(i)]).(leg{l}) = struct;
%         Results.ID.(['session' num2str(i)]).(leg{l}) = struct;
%         Results.SO.(['session' num2str(i)]).(leg{l}) = struct;
%         Results.SO_Activation.(['session' num2str(i)]).(leg{l}) = struct;
%         Results.JRL.(['session' num2str(i)]).(leg{l}) = struct;
%     end
% end


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
% function save_log()


% --------------------------------------------------------------------------------------------------------------- %
function [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct)
% data_struct = a struct similar to the output of load_sto_file



% add mass to results_struct AND remove "


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

        lim1 = round(length(Results.JRL.(session).(leg).HCF)*0.08);
        lim2 = round(length(Results.JRL.(session).(leg).HCF)*0.35);
        lim3 = round(length(Results.JRL.(session).(leg).HCF)*0.8);

        [val1(1,:), loc1(1,:)] = calc_max(Results.JRL.(session).(leg).HCF(lim1:lim2,:),0);
        [val2(1,:), loc2(1,:)] = calc_max(Results.JRL.(session).(leg).HCF(lim2:lim3,:),0);
        loc1(1,:) = loc1(1,:) + lim1-1;
        loc2(1,:) = loc2(1,:) + lim2-1;

        ind = find(val2>val1);
        for k = 1:size(ind,2)
            temp = val1(ind(k));
            val1(ind(k)) = val2(ind(k));
            val2(ind(k)) = temp;
            temp = loc1(ind(k));
            loc1(ind(k)) = loc2(ind(k));
            loc2(ind(k)) = temp;
        end

        Results.JRL.(session).(leg).peak_HCF_val(1,:) = val1;
        Results.JRL.(session).(leg).peak_HCF_loc(1,:) = loc1;
        Results.JRL.(session).(leg).peak_HCF_val(2,:) = val2;
        Results.JRL.(session).(leg).peak_HCF_loc(2,:) = loc2;

        clear val1 val2 loc1 loc2

        %         [Results.JRL.(session).(leg).peak_HCF_val(1,:), Results.JRL.(session).(leg).peak_HCF_loc(1,:)] = calc_max(Results.JRL.(session).(leg).HCF(lim1:lim2,:),0);
        %         [Results.JRL.(session).(leg).peak_HCF_val(2,:), Results.JRL.(session).(leg).peak_HCF_loc(2,:)] = calc_max(Results.JRL.(session).(leg).HCF(lim2:lim3,:),0);
        %         Results.JRL.(session).(leg).peak_HCF_loc(1,:) = Results.JRL.(session).(leg).peak_HCF_loc(1,:) + lim1-1;
        %         Results.JRL.(session).(leg).peak_HCF_loc(2,:) = Results.JRL.(session).(leg).peak_HCF_loc(2,:) + lim2-1;

        % find peak KCF
        Results.JRL.(session).(leg).KCF = sum3D(Results.JRL.(session).(leg), ['knee_' leg(1) '_on_tibia_' leg(1) '_in_tibia_' leg(1)]);

        lim1 = round(length(Results.JRL.(session).(leg).KCF)*0.08);
        lim2 = round(length(Results.JRL.(session).(leg).KCF)*0.35);
        lim3 = round(length(Results.JRL.(session).(leg).KCF)*0.8);

        [val1(1,:), loc1(1,:)] = calc_max(Results.JRL.(session).(leg).KCF(lim1:lim2,:),0);
        [val2(1,:), loc2(1,:)] = calc_max(Results.JRL.(session).(leg).KCF(lim2:lim3,:),0);
        loc1(1,:) = loc1(1,:) + lim1-1;
        loc2(1,:) = loc2(1,:) + lim2-1;

        ind = find(val2>val1);
        for k = 1:size(ind,2)
            temp = val1(ind(k));
            val1(ind(k)) = val2(ind(k));
            val2(ind(k)) = temp;
            temp = loc1(ind(k));
            loc1(ind(k)) = loc2(ind(k));
            loc2(ind(k)) = temp;
        end

        Results.JRL.(session).(leg).peak_KCF_val(1,:) = val1;
        Results.JRL.(session).(leg).peak_KCF_loc(1,:) = loc1;
        Results.JRL.(session).(leg).peak_KCF_val(2,:) = val2;
        Results.JRL.(session).(leg).peak_KCF_loc(2,:) = loc2;

        clear val1 val2 loc1 loc2

        %         [Results.JRL.(session).(leg).peak_KCF_val(1,:), Results.JRL.(session).(leg).peak_KCF_loc(1,:)] = calc_max(Results.JRL.(session).(leg).KCF(lim1:lim2,:),0);
        %         [Results.JRL.(session).(leg).peak_KCF_val(2,:), Results.JRL.(session).(leg).peak_KCF_loc(2,:)] = calc_max(Results.JRL.(session).(leg).KCF(lim2:lim3,:),0);
        %         Results.JRL.(session).(leg).peak_KCF_loc(1,:) = Results.JRL.(session).(leg).peak_KCF_loc(1,:) + lim1-1;
        %         Results.JRL.(session).(leg).peak_KCF_loc(2,:) = Results.JRL.(session).(leg).peak_KCF_loc(2,:) + lim2-1;


        Muscles = fields(Results.SO.(session).(leg));
        for m = 1:length(Muscles)
            muscle = Muscles{m};
            if contains(muscle, {'time', 'durationInSeconds', 'reserve', 'ercspn', 'extobl', 'Nintobl'})
                [Results.SO.(session).(leg).(['peak_' muscle '_val']), Results.SO.(session).(leg).(['peak_' muscle '_loc'])] = calc_max(Results.SO.(session).(leg).(muscle),1);
            else
                [Results.SO.(session).(leg).(['peak_' muscle '_val']),Results.SO.(session).(leg).(['peak_' muscle '_loc'])] = calc_max(Results.SO.(session).(leg).(muscle),0);
            end
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
function [val,loc] = calc_max(time_curve, time_t_f)

% [val, loc] = max(movmean(time_curve,3));
time_curve_smth = movmean(time_curve,3);
for k = 1:size(time_curve_smth,2)
    [val{k}, loc{k}] = findpeaks(time_curve_smth(:,k),'SortStr', 'descend');
    if size(val{k},1) >1
        val{k}(2:end) = [];
        loc{k}(2:end) = [];
    end
    if isempty(val{k}) && time_t_f == 0
        disp('one curve does not have a peak')
        val{k} = NaN;
        loc{k} = NaN;
    elseif isempty(val{k}) && time_t_f == 1
        val{k} = NaN;
        loc{k} = NaN;
    end
end
val = cell2mat(val);
loc = cell2mat(loc);

% --------------------------------------------------------------------------------------------------------------- %
function folders = getfolders(directory, contianing_string, IgnoreCase)

folders = dir (directory);
try
    folders(1:2) =[];                                           % delete "../" and "./"
catch
    return
end
folders = folders([folders.isdir]);                         % select only rows that contain folders

if nargin > 1
    if IgnoreCase == 0
        folders = folders(contains({folders.name},contianing_string));
    else
        folders = folders(contains({folders.name},contianing_string,"IgnoreCase",true));
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,n,Alpha,Color, MakerSize)
% %% Description - Basilio Goncalves (2020)
% https://www.researchgate.net/profile/Basilio_Goncalves
%
% create a scatter plot and a line of best fit with a shaded area as CI
%
%INPUT
%   x = x axis values
%   y = y axis values
%   n = polynomial degree (Default = 1)
%   Alpha = alpha level betwen 0 and 1(Default = 0.05)
%   Color = color of the plots (Default = black)
%   MakerSize = size of the individual markers (Default = 10)
%-------------------------------------------------------------------------
%OUTPUT
%   rsquared = cell matrix maximum Force value
%   pvalue = cell matrix maximum Force value
%   p1 = handle of the plot

if  nargin < 3 || isempty(n)
    n= length(x);
end
if nargin < 4 || isempty(Alpha)
    Alpha= 0.05;
end

if  nargin < 5 || isempty(Color)
    Color= 'k';
end

if  nargin < 6 || isempty(MarkerSize)
    MakerSize= 10;
end

% delete nan
IDXnan = any(isnan([x y]),2);
x(IDXnan)=[];
y(IDXnan)=[];


% calculate polynomial and line that best fit
[p,S,mu] = polyfit(x,y,n);
[y_fit,delta] = polyval(p,x,S,mu);

% calculate rsquared and p-value
[c, pvalue,rlo,rup] = corrcoef(x,y);
rsquared = c(1,2)^2;
pvalue = pvalue(1,2);
rlo = rlo(1,2);
rup = rup(1,2);

% calculate t-statistic
t = tinv(1-Alpha/2,S.df);

% calculate confidence interval
uB = y_fit+ t*delta;
lB = y_fit- t*delta;


p = sortrows([x y y_fit uB lB],1);                                  %plot data
x= p(:,1)';
y= p(:,2)';
y_fit= p(:,3)';
uB= p(:,4)';
lB= p(:,5)';

hold on
p1 = plot(x,y,'o','MarkerSize',MakerSize,'MarkerFaceColor',Color,'MarkerEdgeColor','none');

p1(end+1) = plot(x,y_fit,'Color',Color,'LineWidth',2);              % plot mean trendline

X=[x,fliplr(x)];                                                    % create continuous x value row vector for plotting
Y=[lB fliplr(uB)];                                                  % create y values for out and then back (Remove NaN)
f1 = fill(X,Y,'r');
alpha 0.2                                                           % transparency
set(f1,'FaceColor', Color,'EdgeColor','none')

% --------------------------------------------------------------------------------------------------------------- %
function TimeNormalizedData = TimeNorm (Data,fs)

TimeNormalizedData=[];

for col = 1: size (Data,2)

    currentData = Data(:,col);
    currentData(isnan(currentData))=[];
    if length(currentData)<3
        TimeNormalizedData(1:101,col)= NaN;
        continue
    end

    timeTrial = 0:1/fs:size(currentData,1)/fs;
    timeTrial(end)=[];
    Tnorm = timeTrial(end)/101:timeTrial(end)/101:timeTrial(end);

    TimeNormalizedData(1:101,col)= interp1(timeTrial,currentData,Tnorm)';
end



% add overall title for the figure
suptitle(['muscle forces ' Type])

% add legend
lg = legend({'pre' 'sd' 'post' '' 'td' ''});
lg.FontSize = 12;
lg.Position = [0.8 0.12 0.05 0.1];

% make figure nice (backgorund color, font size and type, etc...)
makeMyFigureNice

% save figure
saveas(gcf, [savedir fp 'muscle_forces' leg '.jpeg'])


% --------------------------------------------------------------------------------------------------------------- %
function plotContactForces(Results) % plot all muscle forces and run t-test SPM to compare (pre,post, & TD)

legs = {'left','right'};
line_colors = getColor('viridis',3);

% JRF names containing fx,fy,fz,fresultant
JRL_names = fields(Results.JRL.session1.left);
JRL_names = [JRL_names(endsWith(JRL_names,{'_fx','_fy','fz','_fresultant'}))]; 

% replace "_l" and "_r_" with "_leg_" to plot just the variable for a
% single leg
JRL_names = strrep(JRL_names,'_r_','_leg_');
JRL_names = strrep(JRL_names,'_l_','_leg_');
JRL_names = unique(JRL_names);

% create figure with subplots
n_subplots = length(JRL_names);
[ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);

% loop through all muscles
for iJRL = 1:length(JRL_names)
    indData = {[] [] []};
    MeanForcesAllSessions = [];
    SDForcesAllSessions = [];
    % assign muscle forces for each session to one
    for iSess = 1:3
        session = ['session' num2str(iSess)];
     
        for iLeg = 1:2
            leg = legs{iLeg};

            % select either absolute or normalised
            contact_forces = Results.JRL.(session).(leg);
            current_JRL_variable =  strrep(JRL_names{iJRL},'_leg_',['_' leg(1) '_']);

            % get all muscle segments for each muscle name
            single_varibale_JRL = contact_forces.(current_JRL_variable);

            % if a column has all zeros make it all NaN
            single_varibale_JRL = ZeroToNaN(single_varibale_JRL);

            % get mass per session (do not move outside the loop to avoid
            % deleting first participant more than once)
            mass_current_session = Results.Mass.(session)(~isnan(Results.Mass.(session)));

            % participant not good data
            if iSess == 3 && length(mass_current_session) > 3
                single_varibale_JRL(:,7)=NaN;
                mass_current_session(1) = [];
                warning on
                warning ('deleting data participant in column 7 for session 3')
            end

            % remove columns with just NaNs
            single_varibale_JRL = removeNaNrows(single_varibale_JRL,2);

            % flip the contact forces for the left leg
            if contains(leg(1),'l') && contains(current_JRL_variable(end-1:end),'fz')
                single_varibale_JRL = -single_varibale_JRL;
            end
            
            % normalise to BW
            single_varibale_JRL_BW_normalised = single_varibale_JRL./(mass_current_session*9.81);

            % find dimensions of the data
            [nRows,nCols] = size(single_varibale_JRL_BW_normalised);
            idxCols = [size(indData{iSess},2)+1 : size(indData{iSess},2)+nCols];

            % add individual participant data to a cell
            indData{iSess}(:,idxCols) = single_varibale_JRL_BW_normalised;

        end

        % calculate mean and SD for each group
        MeanForcesAllSessions(:,iSess) = nanmean(indData{iSess},2);
        SDForcesAllSessions(:,iSess) = nanstd(indData{iSess},0,2);
    end

    parts_JRL_variable = split(current_JRL_variable,'_');
    current_JRL_variable_title = [parts_JRL_variable{1} '_' parts_JRL_variable{end}];

    % run and save SPM plots
    group_types = [1,1,2];
    [SPM] = ttest_spm(indData,group_types);
    suptitle([current_JRL_variable_title])

    savedir = [get_main_dir() fp 'SPM_results' fp 'JRL'];
    if ~isfolder(savedir); mkdir(savedir); end
    saveas(gcf, [savedir fp current_JRL_variable_title '.jpeg'])
    close(gcf)

    % plot force-time curves and add SPM lines on the plot
    axes(ha(iJRL)); hold on
    p = plotShadedSD(MeanForcesAllSessions,SDForcesAllSessions,line_colors);

    % change ylim and yticks
%     ylim([0 6])
    
    add_spm_to_plot(SPM)

    % add ylable to first col
    if any(iJRL == FirstCol)
        yl = ylabel('Joint contact forces (BW)');
        yl.Rotation = 0;
        yl.HorizontalAlignment ="right";
    end

    % title
    t = title(current_JRL_variable_title, 'Interpreter','none');
    t.VerticalAlignment = 'top';
end

% add ticks to the plots
tight_subplot_ticks (ha,LastRow,0)

% add legend
lg = legend({'pre' 'sd' 'post' '' 'td' ''});
lg.FontSize = 12;
lg.Position = [0.8 0.12 0.05 0.1];

% make figure nice (backgorund color, font size and type, etc...)
makeMyFigureNice

% save figure
saveas(gcf, [savedir fp 'JointReactionLoads.jpeg'])



% --------------------------------------------------------------------------------------------------------------- %
function p = plotShadedSD(YData,SD,COLOR,Xvalues)
% plot each column on the Mean matrix with each column on the SD matrix

if size(YData,2)~=size(SD,2)
    error('Number of columans in Mean and SD inputs must agree')
end

% [ cMat, cStruct, cNames] = getColorSet(30); % color blind friendly
if nargin<3 || isempty(COLOR)
    cMat = getColor(0,size(YData,2)); % color scheme 2 (Bas)
else
    cMat = COLOR;
end
%line styles
style = {'-','--',':','-.'};
for k = 1:ceil(size(YData,2)/4)
    style = [style style];
end

for ii = 1:size(YData,2)
    hold on
    y1=YData(:,ii)';                             % create main curve
    if nargin <4
        x=1:length(y1);                             % initialize x row vector
    else
        x = Xvalues';
    end

    p(ii) = plot(x,y1,'LineWidth',1);
    set(p(ii),'Color',cMat(ii,:),'LineStyle',style{ii})
    color = p(ii).Color;

    Top= y1+SD(:,ii)';                          % create top of shaed area
    Bottom = y1-SD(:,ii)';                      % create bottom of shaded
    X=[x,fliplr(x)];                            % create continuous x value array for plotting
    Y=[Bottom fliplr(Top)];                     % create y values for out and then back
    f1 = fill(X,Y,color);
    alpha 0.3
    set(f1,'FaceColor', color,'EdgeColor','none')

end

% --------------------------------------------------------------------------------------------------------------- %

function [ha, pos,FirstCol,LastRow,LastCol] = tight_subplot(Nh, Nw, gap, marg_h, marg_w,Size)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos, FirstCol, LastRow] = tight_subplot(Nh, Nw, gap, marg_h, marg_w,Size)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%       Size     4 dimesions of the figure to be plotted.
%                Size == 0 -> make full size figure
%                Size > 0 (eg 0.4) ->  make size figure 0.4 times full size
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplotBG(3,2,[.01 .03],[.1 .01],[.01 .01],[100 200 1200 600])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
%
% Edited by Basilio Goncalves (2021)
% see also: numSubplots

if nargin < 2 || Nw==0
    RemoveSubPlots  = 1;
    LastPlot        = Nh;
    [N,]            = numSubplots(Nh);
    Nh              = N(1);
    Nw              = N(2);
else
    RemoveSubPlots  = 0;
    LastPlot        = Nh*Nw;
end


if nargin<3 || isempty(gap); gap = .04; end
if nargin<4 || isempty(marg_h); marg_h = .04; end
if nargin<5 || isempty(marg_w); marg_w = .04; end
if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end
figure
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh;
% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);

    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);

if nargin<6 || length(Size)==1 && Size == 0 || length(Size)==1 && Size > 1
    set(gcf,'units','normalized','outerposition',[0 0.05 1 0.95])
elseif length(Size)<3 && Size>0
    dimensions = [1*Size 0.95*Size];
    origin = (0 +(1-dimensions))/2;
    set(gcf,'units','normalized','outerposition',[origin dimensions])
elseif length(Size)==4 && ~any(Size>1)
    set(gcf,'units','normalized','outerposition',[Size])
else
    set(gcf, 'Position', Size);
end

FirstCol = 1:Nw:LastPlot; LastRow = Nh*Nw-Nw+1:Nh*Nw; LastCol = Nw:Nw:LastPlot;

if RemoveSubPlots~=0
    delete(ha(LastPlot+1:end));
    ha(LastPlot+1:end) =[];
    LastRow = ceil(LastPlot/Nw);
    LastRow = LastRow*Nw-Nw+1 : LastPlot;
end

% ------------------------------------------------------------------------- %
function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots.
%
% Inputs
% n - the desired number of subplots.
%
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p =
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end
%
%
% Rob Campbell - January 2010


while isprime(n) & n>4,
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end


% ------------------------------------------------------------------------%

function tight_subplot_ticks (ha,xt,yt)
% tight_subplot_ticks adds xticklabels and yticklables to the assigned
% subplots in "ha" (if xt OR yt == 1, add labels to all subplots)
%
% tight_subplot_ticks (ha,xt,yt)
% xt = index of subplots to add xticks to
% yt = index of subplots to add yticks to
%
% see also tight_subplotBG

N = length(ha);
for i = 1:N
    axes(ha(i));
    if any(xt==i) || any(xt==0)
        xticklabels(xticks)
    else
        xticklabels('')
    end

    if any(yt==i) || any(yt==0)
        yticklabels(yticks)
    else
        yticklabels('')
    end
end


% --------------------------------------------------------------------------------------------------------------- %
function plotKinematics(Results,Type)
if nargin < 2
    Type ='Normalised';
end

legs = {'left','right'};
for iLeg = 1:2
    leg = legs{iLeg};
    IK_values = Results.IK.session1.left;
    IK_names = fields(IK_values);
    
    % only leg side muscles
    IK_names = IK_names(endsWith(IK_names,['_' leg(1)]));
    % remove fileds that contain the term "peaks"
    IK_names = IK_names(~contains(IK_names,{'subtalar','mtp'}));                                                   
    line_colors = getColor('viridis',3);

    % remove '_l' or '_r' from the muscle names
%     IK_names_short = {};
%     for iIK = 1:length(IK_names)
%         IK_names_short{end+1,1} = IK_names{iIK}(1:end-2);
%         if any(contains(IK_names_short{end}(end),{'1','2','3'}))                                                % remove numebers 1,2,3
%             IK_names_short{end} = IK_names_short{end}(1:end-1);
%         end
%     end
%     IK_names_short = unique(IK_names_short);

    % create figure with subplots
    n_subplots = length(IK_names);
    [ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);
  
    % loop through all joints
    for iIK = 1:length(IK_names)
        indData = {};
        MeanForcesAllSessions = [];
        SDForcesAllSessions = [];

        % assign kinematics for each session to one 
        for iSess = 1:3
            session = ['session' num2str(iSess)];

            % select either absolute or normalised
%             switch Type
%                 case 'Normalised'
%                     muscle_forces = Results.SO_normalised.(session).(leg);
%                 case 'Absolute'
%                     IK_values = Results.SO.(session).(leg);
%             end

            
            current_DOF = IK_names{iIK};

            % get all muscle segments for each muscle name
            segments = IK_names(contains(IK_names,current_DOF));                                         
            single_IK = IK_values.(segments{1});

            % if a column has all zeros make it all NaN
            single_IK = ZeroToNaN(single_IK,2);
            
            % average the IK for each segment (for each column /trial)
            for iSeg = 2:length(segments)
                single_IK = (single_IK + IK_values.(segments{iSeg}))./2;                       
            end
            
            if iSess == 3
                single_IK(:,7)=NaN;
            end


            % plot mean and SD
            indData{iSess} = removeNaNrows(single_IK);
            MeanForcesAllSessions(:,iSess) = nanmean(single_IK,2);
            SDForcesAllSessions(:,iSess) = nanstd(single_IK,0,2);

        end

        % run and save SPM plots
        [SPM] = ttest2(indData);
        suptitle([current_DOF ' ' leg])
        
        savedir = [get_main_dir() fp 'SPM_results'];
        if ~isfolder(savedir); mkdir(savedir); end
        saveas(gcf, [savedir fp current_muscle ' ' leg '.jpeg'])
        close(gcf)
        
        % plot force-time curves and add SPM lines on the plot
        axes(ha(iIK)); hold on
        p = plotShadedSD(MeanForcesAllSessions,SDForcesAllSessions,line_colors);
        
        add_spm_to_plot(SPM)
       
        % change ylim and yticks
        ylim([0 140])
        

        % add ylable to first col
        if any(iIK == FirstCol)
            yl = ylabel('% Max isom force');
            yl.Rotation = 0;
            yl.HorizontalAlignment ="right";
        end

        % title
        t = title(current_muscle, 'Interpreter','none');
        t.VerticalAlignment = 'top';
    end
    
    % add ticks to the plots
    tight_subplot_ticks (ha,LastRow,FirstCol)                                                                       
    
    % add overall title for the figure 
    suptitle(['muscle forces ' leg])

    % add legend
    lg = legend({'pre' 'sd' 'post' '' 'td' ''});
    lg.FontSize = 12;
    lg.Position = [0.8 0.12 0.05 0.1];

    % make figure nice (backgorund color, font size and type, etc...)
    makeMyFigureNice
    
    % save figure
    saveas(gcf, [savedir fp 'IK' leg '.jpeg'])

end

% ---------------------------------------------------------------------- %
function [cMat,LineStyles,Marker] = getColor (pallet_name,nColors)

if nColors>163
    error('colorBG function only works for n2 =<63')
end

if nargin < 1
    cMat = convertRGB([176, 104, 16; ...
        16, 157, 176; ...
        136, 16, 176;176,...
        16, 109;31, 28, 28]);  % color scheme 2 (Bas)
else

    switch pallet_name

        case 'prism';   cMat = prism;
        case 'parula';  cMat = parula;
        case 'flag';    cMat = flag;
        case 'hsv';     cMat = hsv;
        case 'hot';     cMat = hot;
        case 'cool';    cMat = cool;
        case 'spring';  cMat = spring;
        case 'summer';  cMat = summer;
        case 'autumn';  cMat = autumn;
        case 'winter';  cMat = winter;
        case 'gray';    cMat = gray;
        case 'bone';    cMat = bone;
        case 'copper';  cMat = copper;
        case 'pink';    cMat = pink;
        case 'lines';   cMat = lines;
        case 'jet';     cMat = jet;
        case 'colorcube';cMat = colorcube;
        case 'viridis'; cMat = viridis;
    end
    warning off
    if nargin == 2
        if nColors>163; error('colorBG function only works for n2 =<63'); end
        cMat = cMat(1:length(cMat)/nColors:length(cMat),:);
    else
        cMat = cMat([1:64],:);
    end

    warning on
end

S = sum(cMat,2);
Duplicates = [];
commonLS = {'-' '--' ':' '-.'};
commonMK = {'none' '+' 's' 'd' '^' 'v'};
N = length(commonMK);
LineStyles ={};Marker ={};
for k = 1:size(cMat,1)
    idx = find(S==S(k));
    if length(idx)==1
        LineStyles{k,1} = commonLS{1};
        Marker{k,1} = commonMK{1};
    else
        if length(idx)/N~= ceil(length(idx)/N)
            idx(end+1:ceil(length(idx)/N)*N)=NaN;
        end
        idx = reshape(idx,N,[]);
        for col = 1:size(idx,2)
            for row = 1:size(idx,1)
                if ~isnan(idx(row,col))
                    LineStyles{idx(row,col),1} = commonLS{col};
                    Marker{idx(row,col),1} = commonMK{row};
                end
            end
        end
    end
    idx = reshape(idx,[],1);idx(isnan(idx))=[];
    S(idx) = NaN;
    Duplicates(idx,1)=k;
end

% -------------------------------------------------------------------- %
function [SPM] = ttest2(indData)

n_groups = length(indData);
combinations = nchoosek([1:n_groups],2); 
Alpha = 0.05;%/size(combinations,1);

ha = tight_subplot(1,n_groups,[],[],[],[0.05 0.3 0.9 0.5]);
SPM = struct;
count = 0;
for iComb = combinations'
    count = count + 1;
    Dataset1 = indData{iComb(1)};
    Dataset2 = indData{iComb(2)};
    
    % run SPM tests
    spmi = spm1d.stats.ttest2(Dataset1',Dataset2');
    spmi = spmi.inference(Alpha);
    
    comparison_name = [num2str(iComb(1)) '_VS_' num2str(iComb(2))];

    % SPM plot with p-values and t-values threshold
    axes(ha(count))
    spmi.plot; spmi.plot_p_values;
    spmi.plot_threshold_label;
    title(comparison_name,'Interpreter','none')

    % from the axis, find the index of the siginifant points and it's p-value
    LinesPlot = ha(count).Children; 
    significant_idx=[]; 
    p_value_idx=[];
    for iLine = 1:length(LinesPlot)

        %  find indexes of patches (signficand shaded areas)   
        if contains(class(LinesPlot(iLine)),'Patch'); significant_idx(end+1) = iLine; end      

        %  find indexes of tesxt with 'P = ' 
        if contains(class(LinesPlot(iLine)),'Text') && contains(LinesPlot(iLine).String,'p '); p_value_idx(end+1) = iLine; end      
    end

    % add spmi results to final struct
    SPM.(['comp_' comparison_name]) = struct;
    flds_spmi = fields(spmi);
    for i = 1:length(flds_spmi)
        curr_fld = flds_spmi{i};
        SPM.(['comp_' comparison_name]).(curr_fld) = spmi.(curr_fld);    
    end

    % add the significance vectors (idx and 
    SPM.(['comp_' comparison_name]).sig_idx = significant_idx;
    SPM.(['comp_' comparison_name]).p_idx = p_value_idx;
end
tight_subplot_ticks(ha,0,0)