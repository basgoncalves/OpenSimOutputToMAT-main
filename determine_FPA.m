function determine_FPA
clear all;
close all;
addpath(genpath('C:\Users\Balu\Nextcloud\Documents\MA\Code\MSKmodelling')); 

% surg = '\pre'; %\pre or \post or ''

subj = struct;
subj.P01_pre = {'dynamic03', 'dynamic04', 'dynamic05', 'dynamic06', 'dynamic07', 'dynamic08', 'dynamic09'};
subj.P01_post = {'Dynamic01', 'Dynamic08', 'Dynamic09', 'Dynamic14'};
subj.P02_pre = {'dynamic06', 'dynamic08', 'dynamic09', 'dynamic27'};
subj.P02_post = {'dynamic07', 'dynamic08', 'dynamic10', 'dynamic12'};
subj.P03_pre = {'Dynamic06', 'Dynamic07', 'Dynamic08', 'Dynamic10', 'Dynamic11'};
subj.P03_post = {'Dynamic07', 'Dynamic08', 'Dynamic09', 'Dynamic13', 'Dynamic14', 'Dynamic15'};
subj.P04_pre = {'Dynamic07', 'Dynamic21', 'Dynamic22', 'Dynamic23', 'Dynamic25'};
subj.P04_post = {'Dynamic05', 'Dynamic07','Dynamic09', 'Dynamic10', 'Dynamic17', 'Dynamic21'};
subj.P05_pre = {'Dynamic5', 'Dynamic7', 'Dynamic9', 'Dynamic11', 'Dynamic12', 'Dynamic17'};
subj.P05_post = {'dynamic09', 'dynamic12', 'dynamic15', 'dynamic21', 'dynamic22'};
subj.TD04 = {'3DGAIT_B_W1', '3DGAIT_B_W2', '3DGAIT_B_W8', '3DGAIT_B_W9', '3DGAIT_B_W10', '3DGAIT_B_W12', '3DGAIT_B_W18', '3DGAIT_B_W19', '3DGAIT_B_W20' }; %TD04
subj.TD06 = {'3DGAIT_B_W12', '3DGAIT_B_W19', '3DGAIT_B_W14', '3DGAIT_B_W20'}; %TD06
subj.TD07 = {'3DGAIT_A_W6', '3DGAIT_A_W7', '3DGAIT_A_W9', '3DGAIT_A_W18', '3DGAIT_A_W25', '3DGAIT_A_W27', '3DGAIT_A_W29', '3DGAIT_A_W32'}; %TD07

fields_names = fieldnames(subj);
fpa_av_l = nan(1,length(fields_names));
fpa_av_r = nan(1,length(fields_names));
fpa_char = struct;
%% calculate body kinematics

% edit xml file
% temp = xml_read('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\P01\pre\bodyKinematics_template.xml');
% temp.AnalyzeTool.ControllerSet = [];
% temp.AnalyzeTool.final_time = 'Inf';

for k=1:length(fields(subj))
    if startsWith(fields_names{k},'P')
        subj_id = fields_names{k}(1:3);
        surg = erase(fields_names{k},[subj_id '_']);
        surg = ['\' surg];
    else
        subj_id = fields_names{k}(1:4);
        surg = '';
    end

% temp.AnalyzeTool.model_file = ['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\Model\FINAL_PERSONALISEDTORSIONS_scaled_final.osim'];
% 
% % loop through trials
% for i=1:size(subj.(fields_names{k}),2)
%     temp.AnalyzeTool.results_directory = ['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\' subj.(fields_names{k}){i} '\Output\BodyKin'];
%     temp.AnalyzeTool.coordinates_file = ['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\output automization/FINAL_PERSONALISEDTORSIONS_scaled_final/' subj.(fields_names{k}){i} '/Output/IK/IK.mot'];
%     temp.AnalyzeTool.AnalysisSet.objects.BodyKinematics.on = 'true';
%     temp.AnalyzeTool.AnalysisSet.objects.BodyKinematics.in_degrees = 'true';
%     temp.AnalyzeTool.AnalysisSet.objects.BodyKinematics.express_results_in_body_local_frame = 'false';
%     % save xml file
%     root = 'OpenSimDocument';
%     Pref = struct;
%     Pref.StructItem = false;
%     Pref.CellItem = false;
%     temp = ConvertLogicToString(temp);
%     path_xml_file = ['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\output automization/FINAL_PERSONALISEDTORSIONS_scaled_final/' subj.(fields_names{k}){i} '/Output/BodyKin/BodyKinematics_setup.xml'];
%     xml_write(path_xml_file, temp, root,Pref);
% 
%     % run analysis using Opensim API
%     import org.opensim.modeling.*
%     osimModel   = Model(temp.AnalyzeTool.model_file);
%     analyzeTool = AnalyzeTool(path_xml_file);
%     analyzeTool.setResultsDir(temp.AnalyzeTool.results_directory);
%     analyzeTool.setModel(osimModel);
%     analyzeTool.setName(' ')
%     analyzeTool.run;
% end

%% calculate FPA
fpa = struct;
fpa.left = [];
l=1;
fpa.right = [];
r=1;
load(['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\output automization\bodyKinematics.mat']);
fields_bodyKin = fields(data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final);

for i=1:size(fields_bodyKin,1) % loop through steps
    %     load(['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\' subj.(fields_names{k}){i} '\settings.mat'])
    %     bodyKin = load_sto_file(['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\' subj_id surg '\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\' subj.(fields_names{k}){i} '\Output\BodyKin\_BodyKinematics_pos_global.sto']);

    if contains(fields_bodyKin{i},'left')
        %         for j=1:size(cycle.left.start,2)
        calc_stance = data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final.(fields_bodyKin{i}).calcn_l_Oy;
        %         calc_stance = bodyKin.calcn_l_Oy(cycle.left.start(j)+1:cycle.left.footOff(j));
        fpa.left(l) = mean(calc_stance);
        l=l+1;
        %         end
    end
    if contains(fields_bodyKin{i},'right')
        %         for j=1:size(cycle.right.start,2)
        calc_stance = data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final.(fields_bodyKin{i}).calcn_r_Oy;
        %         calc_stance = bodyKin.calcn_r_Oy(cycle.right.start(j)+1:cycle.right.footOff(j));
        fpa.right(r) = -mean(calc_stance);
        r=r+1;
        %         end
    end
end

% mean for every participant
if ~isempty(fpa.left)
    fpa_char.left(k) = mean(fpa.left);
end
if ~isempty(fpa.right)
    fpa_char.right(k) = mean(fpa.right);
end


end

% fpa_char.left(11:13) = - fpa_char.left(11:13);
% fpa_char.right(11:13) = - fpa_char.right(11:13);
%% plot results
figure
plot(fpa_char.left(1:5),'*k')
hold on
plot(fpa_char.right(1:5),'*k')
plot(fpa_char.right(6:10),'*r')
plot(fpa_char.left(6:10),'*r')
plot(fpa_char.left(11:13),'*b')
plot(fpa_char.right(11:13),'*b')
legend({'pre', '', 'post', '', 'TD'})