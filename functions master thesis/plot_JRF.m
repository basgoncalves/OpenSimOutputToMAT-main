function plot_JRF(trial_list, joints_var_list, subject_folder, events_folder, file_path_end, y_label, subj_mass, BW_norm, f)
[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(length(joints_var_list),0);

x_val_norm_l = 0;
x_val_norm_r = 0;
nb_l = 0;
nb_r = 0;
m = 1;
n = 1;

for i = 1:length(trial_list)
    data = load_sto_file([subject_folder trial_list{i} file_path_end]);

    ev = importdata([events_folder trial_list{i} '\settings.mat']);
    if isfield(ev.cycle,'left')
        for k = 1:length(ev.cycle.left.start)
            L_GC{k} = ev.cycle.left.start(k);
            if L_GC{k} == 0
                L_GC{k} = 1;
            end
            L_GC{k}(2) = ev.cycle.left.end(k);
            nb_l = nb_l +1;
        end 
    end
    if isfield(ev.cycle,'right')
        for k = 1:length(ev.cycle.right.start)
            R_GC{k} = ev.cycle.right.start(k);
            if R_GC{k} == 0
                R_GC{k} = 1;
            end
            R_GC{k}(2) = ev.cycle.right.end(k);
            nb_r = nb_r +1;
        end
    end

    fs = 1 / (data.time(2)-data.time(1));

    % plot right leg
    if isfield(ev.cycle,'right')
        for iDOF = 1:length(joints_var_list)/2
            axes(ax(iDOF)); hold on; grid on
            name_data_fx = [joints_var_list{1,iDOF} '_fx'];
            name_data_fy = [joints_var_list{1,iDOF} '_fy'];
            name_data_fz = [joints_var_list{1,iDOF} '_fz'];
            data_fx = data.(name_data_fx)(R_GC{1}(1):R_GC{1}(2));
            data_fy = data.(name_data_fy)(R_GC{1}(1):R_GC{1}(2));
            data_fz = data.(name_data_fz)(R_GC{1}(1):R_GC{1}(2));
            raw_data_r = sqrt(data_fx.^2 + data_fy.^2 + data_fz.^2); % resultant force

%             raw_data_r = data.(joints_var_list{iDOF})(R_GC{1}(1):R_GC{1}(2));
            time_norm_data_r = TimeNorm(raw_data_r,fs);
            if length(ev.cycle.right.start) > 1
                for k = 2:length(ev.cycle.right.start)
                    data_fx = data.(name_data_fx)(R_GC{k}(1):R_GC{k}(2));
                    data_fy = data.(name_data_fy)(R_GC{k}(1):R_GC{k}(2));
                    data_fz = data.(name_data_fz)(R_GC{k}(1):R_GC{k}(2));
                    raw_data_r = sqrt(data_fx.^2 + data_fy.^2 + data_fz.^2);
                    raw_data_r = TimeNorm(raw_data_r,fs);
                    time_norm_data_r = time_norm_data_r + raw_data_r;
                end
                time_norm_data_r = time_norm_data_r./length(ev.cycle.right.start);
            end
            if BW_norm == 1
                weight_norm_data_r = time_norm_data_r./(subj_mass*9.81);
                plot(weight_norm_data_r);
%                 data_r{iDOF}(m,:) = weight_norm_data_r;
            else
                plot(time_norm_data_r);
%                 data_r{iDOF}(m,:) = time_norm_data_r;
            end
            title(joints_var_list{iDOF}, 'Interpreter', 'none')
%             if i == length(trial_list)
%             legend(trial_list{1,:}, 'toe off', 'Interpreter', 'None')
%             end
            if any(iDOF == FirstCol)
                ylabel([y_label ' of the right leg'])
            end
        end
        m=m+1;
        %add line indicating right toe off (average of all trials)
        x_val = ev.cycle.right.footOff-ev.cycle.right.start;
        x_val_norm_r = x_val_norm_r + (x_val / (ev.cycle.right.end - ev.cycle.right.start)*100); % time normalize
    end

    % plot left leg
    if isfield(ev.cycle,'left')
        for iDOF = length(joints_var_list)/2+1:length(joints_var_list)
            axes(ax(iDOF)); hold on; grid on
            name_data_fx = [joints_var_list{1,iDOF} '_fx'];
            name_data_fy = [joints_var_list{1,iDOF} '_fy'];
            name_data_fz = [joints_var_list{1,iDOF} '_fz'];
            data_fx = data.(name_data_fx)(L_GC{1}(1):L_GC{1}(2));
            data_fy = data.(name_data_fy)(L_GC{1}(1):L_GC{1}(2));
            data_fz = data.(name_data_fz)(L_GC{1}(1):L_GC{1}(2));
            raw_data_l = sqrt(data_fx.^2 + data_fy.^2 + data_fz.^2);

%             raw_data_r = data.(joints_var_list{iDOF})(R_GC{1}(1):R_GC{1}(2));
            time_norm_data_l = TimeNorm(raw_data_l,fs);
            if length(ev.cycle.left.start) > 1
                for k = 2:length(ev.cycle.left.start)
                    data_fx = data.(name_data_fx)(L_GC{k}(1):L_GC{k}(2));
                    data_fy = data.(name_data_fy)(L_GC{k}(1):L_GC{k}(2));
                    data_fz = data.(name_data_fz)(L_GC{k}(1):L_GC{k}(2));
                    raw_data_l = sqrt(data_fx.^2 + data_fy.^2 + data_fz.^2);
                    raw_data_l = TimeNorm(raw_data_l,fs);
                    time_norm_data_l = time_norm_data_l + raw_data_l;
                end
                time_norm_data_l = time_norm_data_l./length(ev.cycle.left.start);
            end
            if BW_norm == 1
                weight_norm_data_l = time_norm_data_l./(subj_mass*9.81);
%                 data_l{iDOF-length(joints_var_list)/2}(n,:)  = weight_norm_data_l;
                plot(weight_norm_data_l);
            else
%                 data_l{iDOF-length(joints_var_list)/2}(n,:)  = time_norm_data_l;
                plot(time_norm_data_l);
            end
            title(joints_var_list{iDOF}, 'Interpreter', 'none')
%             if i == length(trial_list)
%             legend(trial_list{1,:}, 'toe off', 'Interpreter', 'None')
%             end
            if any(iDOF == FirstCol)
                ylabel([y_label ' of the left leg'])
            end
        end
        n=n+1;
        %add line indicating left toe off (average of all trials)
        x_val = ev.cycle.left.footOff-ev.cycle.left.start;
        x_val_norm_l = x_val_norm_l + (x_val / (ev.cycle.left.end - ev.cycle.left.start)*100); % time normalize
    end

%     for i = 1: length(joints_var_list)
%         axes(ax(i)); hold on; grid on
%         name_data_fx = [joints_var_list{1,i} '_fx'];
%         name_data_fy = [joints_var_list{1,i} '_fy'];
%         name_data_fz = [joints_var_list{1,i}{1} '_fz'];
%         data_fx = data.(name_data_fx)(hs1:hs2,:)/(subj_mass*9.81);
%         data_fy = data.(name_data_fy)(hs1:hs2,:)/(subj_mass*9.81);
%         data_fz = data.(name_data_fz)(hs1:hs2,:)/(subj_mass*9.81);
% 
%         data_res = sqrt(data_fx.^2 + data_fy.^2 + data_fz.^2);
%         data_time_norm = TimeNorm (data_res, fs);
%         plot(data_time_norm)
%         plotVert(to_val_norm)
%         title(var_list{2,i}, 'Interpreter', 'none')
%         if any(i == FirstCol)
%             ylabel(['resultant JRF of the ' s_name ' leg [N/BW]'])
%         end
%     end
% end

end

% add vertical line indication toe off
if isfield(ev.cycle,'right')
    x_val_norm_r = x_val_norm_r/nb_r;
    for k = 1:length(ax)/2
        axes(ax(k));
        plotVert(x_val_norm_r);
    end
end

if isfield(ev.cycle,'left')
    x_val_norm_l = x_val_norm_l/nb_r;
    for k = length(ax)/2+1:length(ax)
        axes(ax(k));
        plotVert(x_val_norm_l);
    end
end
tight_subplot_ticks(ax,LastRow,0)
mmfn

end
% end