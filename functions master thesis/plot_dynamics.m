clear all; close all; clc
% import org.opensim.modeling.*
subj_mass = 32.9;

trial_list1 = {'dyn_03', 'dyn_04', 'dyn_05', 'dyn_06', 'dyn_07', 'dyn_08', 'dyn_09'};
trial_list2 = {'dynamic03'};%, 'dynamic04', 'dynamic05', 'dynamic06', 'dynamic07', 'dynamic08', 'dynamic09'};
% % angles_list = {'hip flexion right', 'hip flexion left', 'hip adduction right', 'hip adduction left', 'hip rotation right', 'hip rotation left', 'knee angle right', 'knee angle left', 'ankle angle right', 'ankle angle left'};
angles_var_list = {'hip_flexion_r_moment',  'knee_angle_r_moment',  'ankle_angle_r_moment', 'hip_flexion_l_moment','knee_angle_l_moment','ankle_angle_l_moment'};

subject_folder = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\C3D\';

[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(6,0);

for i = 1:length(trial_list2)

%     dyn = load_sto_file([subject_folder trial_list2{i} '\inverse_dynamics.sto']);
    dyn = load_sto_file('C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\dynamic03\Output\ID\inverse_dynamics.sto');
    ev = importdata([subject_folder trial_list2{i} '\settings.mat']);
    if isfield(ev.cycle,'left')
        L_GC = ev.cycle.left.start(1);
        L_GC(2) = ev.cycle.left.end(1);
    end
    if isfield(ev.cycle,'right')
        R_GC = ev.cycle.right.start(1);
        R_GC(2) = ev.cycle.right.end(1);
    end

    fs = 1 / (dyn.time(2)-dyn.time(1));
 
    % plot right leg
    if isfield(ev.cycle,'right')
        for iDOF = 1:length(angles_var_list)/2
            axes(ax(iDOF)); hold on; grid on
    
            raw_data = dyn.(angles_var_list{iDOF})(R_GC(1):R_GC(2));
            time_norm_data = TimeNorm(raw_data,fs);
            weight_norm_data = time_norm_data/(subj_mass);
            plot(weight_norm_data);
            title(angles_var_list{iDOF}, 'Interpreter', 'none')
            if any(iDOF == FirstCol)
                ylabel('Joint moments (Nm)')
            end
        end
    end
    % plot left leg
    if isfield(ev.cycle,'left')
        for iDOF = length(angles_var_list)/2+1:length(angles_var_list)
            axes(ax(iDOF)); hold on; grid on
    
            raw_data = dyn.(angles_var_list{iDOF})(L_GC(1):L_GC(2));
            time_norm_data = TimeNorm(raw_data,fs);
            weight_norm_data = time_norm_data/(subj_mass);
            plot(weight_norm_data);
            title(angles_var_list{iDOF}, 'Interpreter', 'none')
            if any(iDOF == FirstCol)
                ylabel('Joint moments (Nm)')
            end
        end
    end
end

%add line indicating right toe off (of last available trial)
x_val = ev.cycle.right.footOff-ev.cycle.right.start;
x_val_norm = x_val / (ev.cycle.right.end - ev.cycle.right.start)*100; % time normalize
for k = 1:length(ax)/2
    axes(ax(k));
    plotVert(x_val_norm);
end
% add line indicating left toe off (of last available trial)
x_val = ev.cycle.left.footOff-ev.cycle.left.start;
x_val_norm = x_val / (ev.cycle.left.end - ev.cycle.left.start)*100; % time normalize
for k = length(ax)/2+1:length(ax)
    axes(ax(k));
    plotVert(x_val_norm);
end

tight_subplot_ticks (ax,LastRow,0)

mmfn