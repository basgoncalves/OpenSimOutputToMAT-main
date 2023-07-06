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

% data_pre = struct;
% data_post = struct;
% data_TD = struct;

load('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\Results_BW.mat');
muscles_r = append(muscle_groups,'_r');
muscles_l = append(muscle_groups,'_l');

alpha_jcf = alpha;
alpha_m = alpha;

%% JCF

% for i = 1:size(joints,1)
%     data_pre = [Results_BW.JRL.(sessions{1}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{1}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
%     data_post = [Results_BW.JRL.(sessions{2}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{2}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
%     data_TD = [Results_BW.JRL.(sessions{3}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{3}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
% 
% % paired: data_pre vs data_post
% diff_pre_post = data_pre - data_post; % neg if values higher post surg; positive when values are lower post surg
% [jcf.paired_t.h(i),jcf.paired_t.p(i),jcf.paired_t.ci(i,:),jcf.paired_t.stats(i)] = ttest(diff_pre_post,0,alpha,'both');
% 
% % two sample: data_pre vs data_td && data_post vs data_td
% [jcf.two_t_pre.h(i),jcf.two_t_pre.p(i),jcf.two_t_pre.ci(i,:),jcf.two_t_pre.stats(i)] = ttest2(data_TD, data_pre,alpha); % '2' because two sample t-test
% [jcf.two_t_post.h(i),jcf.two_t_post.p(i),jcf.two_t_post.ci(i,:),jcf.two_t_post.stats(i)] = ttest2(data_TD, data_post,alpha); % '2' because two sample t-test
% end

%% muscle forces
% for i = 1:size(muscles_l,1)
%     data_pre = [Results_BW.SO.(sessions{1}).(legs{1}).(muscles_l{i})(1,:), Results_BW.SO.(sessions{1}).(legs{2}).(muscles_r{i})(1,:)];
%     data_post = [Results_BW.SO.(sessions{2}).(legs{1}).(muscles_l{i})(1,:), Results_BW.SO.(sessions{2}).(legs{2}).(muscles_r{i})(1,:)];
%     data_TD = [Results_BW.SO.(sessions{3}).(legs{1}).(muscles_l{i})(1,:), Results_BW.SO.(sessions{3}).(legs{2}).(muscles_r{i})(1,:)];
% 
% % paired: data_pre vs data_post
% diff_pre_post = data_pre - data_post; % neg if values higher post surg; positive when values are lower post surg
% [m.paired_t.h(i),m.paired_t.p(i),m.paired_t.ci(i,:),m.paired_t.stats(i)] = ttest(diff_pre_post,0,alpha_m,'both');
% 
% % two sample: data_pre vs data_td && data_post vs data_td
% [m.two_t_pre.h(i),m.two_t_pre.p(i),m.two_t_pre.ci(i,:),m.two_t_pre.stats(i)] = ttest2(data_TD, data_pre,alpha_m); % '2' because two sample t-test
% [m.two_t_post.h(i),m.two_t_post.p(i),m.two_t_post.ci(i,:),m.two_t_post.stats(i)] = ttest2(data_TD, data_post,alpha_m); % '2' because two sample t-test
% end
% disp('')

%% walking velocity pre vs post
% vel_pre = [1.37, 1.25, 1.38, 1.28, 1.03];
% vel_post = [1.19, 1.21, 1.42, 1.15, 0.97];
% 
% diff_pre_post = vel_pre - vel_post; % neg if values higher post surg; positive when values are lower post surg
% [vel.paired_t.h, vel.paired_t.p, vel.paired_t.ci(:), vel.paired_t.stats] = ttest(diff_pre_post,0,alpha,'both');

%% deviation in joint loads pre vs post
for i = 1:size(joints,1)
    data_pre.(joints{i}) = [Results_BW.JRL.(sessions{1}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{1}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
    data_post.(joints{i}) = [Results_BW.JRL.(sessions{2}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{2}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
    data_TD.(joints{i}) = [Results_BW.JRL.(sessions{3}).(legs{1}).(['peak_' joints{i} '_val'])(1,:), Results_BW.JRL.(sessions{3}).(legs{2}).(['peak_' joints{i} '_val'])(1,:)];
    av_TD.(joints{i}) = mean(data_TD.(joints{i}));
    diff_pre.(joints{i}) = data_pre.(joints{i}) - av_TD.(joints{i});
    diff_post.(joints{i}) = data_post.(joints{i}) - av_TD.(joints{i});
    diff_pre_mean(i) = mean(diff_pre.(joints{i}));
    diff_post_mean(i) = mean(diff_post.(joints{i}));
    diff_pre_sd(i) = std(diff_pre.(joints{i}));
    diff_post_sd(i) = std(diff_post.(joints{i}));

    % paired: diff_pre vs diff_post
    diff_pre_post = diff_pre.(joints{i}) - diff_post.(joints{i}); % neg if values higher post surg; positive when values are lower post surg
    [dev_jcf.paired_t.h(i),dev_jcf.paired_t.p(i),dev_jcf.paired_t.ci(i,:),dev_jcf.paired_t.stats(i)] = ttest(diff_pre_post,0,alpha,'both');
end

%% FPA


% -------------------- FUNCTIONS --------------------------------------%

function data_path = get_main_dir()
activeFile = [mfilename('fullpath') '.m'];
data_path = [fileparts(fileparts(fileparts(activeFile))) '\Kira_MSc_data'];
