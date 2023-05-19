function plot_dyn(trial_list, moments_var_list, subject_folder, events_folder, f)

[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(length(moments_var_list),0);

x_val_norm_l = 0;
x_val_norm_r = 0;

for i = 1:length(trial_list)

    dyn = load_sto_file([subject_folder trial_list{i} '\Output\ID\inverse_dynamics.sto']);
    ev = importdata([events_folder trial_list{i} '\settings.mat']);
    if isfield(ev.cycle,'left')
        for k = 1:length(ev.cycle.left.start)
            L_GC{k} = ev.cycle.left.start(k);
            L_GC{k}(2) = ev.cycle.left.end(k);
        end  
    end
    if isfield(ev.cycle,'right')
        for k = 1:length(ev.cycle.right.start)
            R_GC{k} = ev.cycle.right.start(k);
            R_GC{k}(2) = ev.cycle.right.end(k);
        end
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
end