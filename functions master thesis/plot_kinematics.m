clear all; close all; clc;
% import org.opensim.modeling.*


trial_list1 = {'dyn_03', 'dyn_04', 'dyn_05', 'dyn_06', 'dyn_07', 'dyn_08', 'dyn_09'};
trial_list2 = {'dynamic03', 'dynamic04', 'dynamic05', 'dynamic06', 'dynamic07', 'dynamic08', 'dynamic09'};
% angles_list = {'hip flexion right', 'hip flexion left', 'hip adduction right', 'hip adduction left', 'hip rotation right', 'hip rotation left', 'knee angle right', 'knee angle left', 'ankle angle right', 'ankle angle left'};
angles_var_list = {'hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r',  'knee_angle_r',...
    'ankle_angle_r', 'hip_flexion_l','hip_adduction_l','hip_rotation_l','knee_angle_l','ankle_angle_l'};

subject_folder = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\C3D\';

[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(10,0);

for i = 1:length(trial_list2)


    kin = load_sto_file([subject_folder trial_list2{i} '\IK_Results.mot']);
    ev = importdata([subject_folder trial_list2{i} '\settings.mat']);
    if isfield(ev.cycle,'left')
        L_GC = ev.cycle.left.start(1);
        L_GC(2) = ev.cycle.left.end(1);
    end
    if isfield(ev.cycle,'right')
        R_GC = ev.cycle.right.start(1);
        R_GC(2) = ev.cycle.right.end(1);
    end

    fs = 1 / (kin.time(2)-kin.time(1));
 
    % plot right leg
    if isfield(ev.cycle,'right')
        for iDOF = 1:length(angles_var_list)/2
            axes(ax(iDOF)); hold on; grid on
    
            raw_data = kin.(angles_var_list{iDOF})(R_GC(1):R_GC(2));
            time_norm_data = TimeNorm(raw_data,fs);
            plot(time_norm_data);
            title(angles_var_list{iDOF}, 'Interpreter', 'none')
            if any(iDOF == FirstCol)
                ylabel('Joint angles (deg)')
            end
        end
    end
    % plot left leg
    if isfield(ev.cycle,'left')
        for iDOF = length(angles_var_list)/2+1:length(angles_var_list)
            axes(ax(iDOF)); hold on; grid on
    
            raw_data = kin.(angles_var_list{iDOF})(L_GC(1):L_GC(2));
            time_norm_data = TimeNorm(raw_data,fs);
            plot(time_norm_data);
            title(angles_var_list{iDOF}, 'Interpreter', 'none')
            if any(iDOF == FirstCol)
                ylabel('Joint angles (deg)')
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