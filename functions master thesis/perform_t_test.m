% ------------------------------------------------------------------------%
function [paired_t, two_t_pre, two_t_post]=perform_t_test(calc)
%h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise
%confidence interval ci for the mean of x, or of x – y for the paired t-test, and the structure stats containing information about the test statistic.
alpha = 0.05;
m = struct;
jcf = struct;
% paired_t = struct;
% two_t_pre = struct;
% two_t_post = struct;
legs = {'left'; 'right'};
sessions = {'session1'; 'session2'; 'session3'};
joints = {'HCF'; 'KCF'; 'ACF'};
muscle_groups = {'glut_max'; 'glut_med'; 'glut_min';'adductors'; 'hamstrings'; 'iliopsoas'; 'tfl'; 'ext_rot'; 'rect_fem'};
muscles_r = append(muscle_groups,'_r');
muscles_l = append(muscle_groups,'_l');

data = struct;
% data_post = struct;
% data_TD = struct;

load('C:\Users\Balu\Nextcloud\Documents\MA\Code\OpenSimOutputToMAT-main\Results_BW.mat');
muscles_r = append(muscle_groups,'_r');
muscles_l = append(muscle_groups,'_l');

alpha_jcf = alpha;
alpha_m = alpha;

%% JCF

for i = 1:size(joints,1)
    data_pre = [Results_BW.JRL.(sessions{1}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{1}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
    data_post = [Results_BW.JRL.(sessions{2}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{2}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
    data_TD = [Results_BW.JRL.(sessions{3}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{3}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
data.(joints{i}).mean(1) = mean(data_pre);
data.(joints{i}).mean(2) = mean(data_post);
data.(joints{i}).mean(3) = mean(data_TD);
data.(joints{i}).SD(1) = std(data_pre);
data.(joints{i}).SD(2) = std(data_post);
data.(joints{i}).SD(3) = std(data_TD);


% paired: data_pre vs data_post
diff_pre_post = data_pre - data_post; % neg if values higher post surg; positive when values are lower post surg
[jcf.paired_t.h(i),jcf.paired_t.p(i),jcf.paired_t.ci(i,:),jcf.paired_t.stats(i)] = ttest(diff_pre_post,0,alpha_jcf,'both');

% two sample: data_pre vs data_td && data_post vs data_td
[jcf.two_t_pre.h(i),jcf.two_t_pre.p(i),jcf.two_t_pre.ci(i,:),jcf.two_t_pre.stats(i)] = ttest2(data_TD, data_pre,alpha_jcf); % '2' because two sample t-test
[jcf.two_t_post.h(i),jcf.two_t_post.p(i),jcf.two_t_post.ci(i,:),jcf.two_t_post.stats(i)] = ttest2(data_TD, data_post,alpha_jcf); % '2' because two sample t-test
end
clear data_pre data_post data_TD
%% muscle forces
for i = 1:size(muscles_l,1)
    data_pre = [Results_BW.SO.(sessions{1}).(legs{1}).(['peak_' muscles_l{i} '_val'])(1,:), Results_BW.SO.(sessions{1}).(legs{2}).(['peak_' muscles_r{i} '_val'])(1,:)];
    data_post = [Results_BW.SO.(sessions{2}).(legs{1}).(['peak_' muscles_l{i} '_val'])(1,:), Results_BW.SO.(sessions{2}).(legs{2}).(['peak_' muscles_r{i} '_val'])(1,:)];
    data_TD = [Results_BW.SO.(sessions{3}).(legs{1}).(['peak_' muscles_l{i} '_val'])(1,:), Results_BW.SO.(sessions{3}).(legs{2}).(['peak_' muscles_r{i} '_val'])(1,:)];

    data.(muscle_groups{i}).mean(1) = mean(data_pre);
    data.(muscle_groups{i}).mean(2) = mean(data_post);
    data.(muscle_groups{i}).mean(3) = mean(data_TD);
    data.(muscle_groups{i}).SD(1) = std(data_pre);
    data.(muscle_groups{i}).SD(2) = std(data_post);
    data.(muscle_groups{i}).SD(3) = std(data_TD);

% paired: data_pre vs data_post
diff_pre_post = data_pre - data_post; % neg if values higher post surg; positive when values are lower post surg
[m.paired_t.h(i),m.paired_t.p(i),m.paired_t.ci(i,:),m.paired_t.stats(i)] = ttest(diff_pre_post,0,alpha_m,'both');

% two sample: data_pre vs data_td && data_post vs data_td
[m.two_t_pre.h(i),m.two_t_pre.p(i),m.two_t_pre.ci(i,:),m.two_t_pre.stats(i)] = ttest2(data_TD, data_pre,alpha_m); % '2' because two sample t-test
[m.two_t_post.h(i),m.two_t_post.p(i),m.two_t_post.ci(i,:),m.two_t_post.stats(i)] = ttest2(data_TD, data_post,alpha_m); % '2' because two sample t-test
end

save('data_t_test_peaks.mat', 'data')
save('results_t_test_muscles.mat', 'm')
save('results_t_test_jcf.mat', 'jcf')

disp('')

% -------------------- FUNCTIONS --------------------------------------%

function data_path = get_main_dir()
activeFile = [mfilename('fullpath') '.m'];
data_path = [fileparts(fileparts(fileparts(activeFile))) '\Kira_MSc_data'];
