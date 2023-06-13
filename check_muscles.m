data = load('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\P01\post\output automization\dataStruct_ErrorScores_no_trials_removed.mat');
f = fields(data.data.SO.FINAL_PERSONALISEDTORSIONS_scaled_final);
f_trials = f(~contains(f,'mass'));
f_muscles = fields(data.data.SO.FINAL_PERSONALISEDTORSIONS_scaled_final.(f_trials{1}));
f_muscles = f_muscles(~contains(f_muscles,{'time', 'reserve','duration','FX','FY', 'FZ', 'MX', 'MY', 'MZ', 'calcn'}));
f_muscles_l = f_muscles(endsWith(f_muscles,'_l'));
f_muscles_r = f_muscles(endsWith(f_muscles,'_r'));


for i = 1:length(f_trials)
    figure()
    hold on
    title(f_trials{i})
    [ax, pos,FirstCol,LastRow,LastCol] = tight_subplotBG(length(f_muscles),0);
    for l = 1:length(f_muscles_l)
        axes(ax(l)); hold on; grid on;
        if contains(f_trials{i},'left')
        plot(data.data.SO.FINAL_PERSONALISEDTORSIONS_scaled_final.(f_trials{i}).(f_muscles_l{l}))
        title(f_muscles_l{l})
        else
plot(data.data.SO.FINAL_PERSONALISEDTORSIONS_scaled_final.(f_trials{i}).(f_muscles_r{l}))
        title(f_muscles_r{l})
        end
    end
end