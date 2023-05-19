function [data_r, data_l] = collect_data(trial_list, var_list, subject_folder, events_folder, file_path_end, y_label, subj_mass, BW_norm, f)
m = 1;
n = 1;
data_r = [];
data_l = [];
glut_max = 0;


% [ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(length(var_list),0);

x_val_norm_l = 0;
x_val_norm_r = 0;
nb_l = 0;
nb_r = 0;

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
        for iDOF = 1:length(var_list)/2
%             axes(ax(iDOF)); hold on; grid on
%             switch var_list{iDOF}
%                 case 'glut_max'
%                     glut_max = 1;
%             end
%             if glut_max == 0
                raw_data_r = data.(var_list{iDOF})(R_GC{1}(1):R_GC{1}(2));
                time_norm_data_r = TimeNorm(raw_data_r,fs);
                if length(ev.cycle.right.start) > 1
                    for k = 2:length(ev.cycle.right.start)
                        raw_data_r = TimeNorm(data.(var_list{iDOF})(R_GC{k}(1):R_GC{k}(2)),fs);
                        time_norm_data_r = time_norm_data_r + raw_data_r;
                    end
                    time_norm_data_r = time_norm_data_r./length(ev.cycle.right.start);
                end
%             elseif glut_max == 1
%                 raw_data1_r = data.([var_list{iDOF} '1_r'])(R_GC{1}(1):R_GC{1}(2));
%                 raw_data2_r = data.([var_list{iDOF} '2_r'])(R_GC{1}(1):R_GC{1}(2));
%                 raw_data3_r = data.([var_list{iDOF} '3_r'])(R_GC{1}(1):R_GC{1}(2));
%                 raw_data_r = sqrt(raw_data1_r.^2 + raw_data2_r.^2 + raw_data3_r.^2);
%                 time_norm_data_r = TimeNorm(raw_data_r,fs);
%                 if length(ev.cycle.right.start) > 1
%                     for k = 2:length(ev.cycle.right.start)
%                         raw_data1_r = data.([var_list{iDOF} '1_r'])(R_GC{1}(1):R_GC{1}(2));
%                         raw_data2_r = data.([var_list{iDOF} '2_r'])(R_GC{1}(1):R_GC{1}(2));
%                         raw_data3_r = data.([var_list{iDOF} '3_r'])(R_GC{1}(1):R_GC{1}(2));
%                         raw_data_r = sqrt(raw_data1_r.^2 + raw_data2_r.^2 + raw_data3_r.^2);
%                         raw_data_r = TimeNorm(raw_data_r,fs);
%                         time_norm_data_r = time_norm_data_r + raw_data_r;
%                     end
%                     time_norm_data_r = time_norm_data_r./length(ev.cycle.right.start);
%                 end
%             end
            if BW_norm == 1
                weight_norm_data_r = time_norm_data_r./(subj_mass*9.81);
                data_r{iDOF}(m,:) = weight_norm_data_r;
                %plot(weight_norm_data_r);
            else
                data_r{iDOF}(m,:) = time_norm_data_r;
                %plot(time_norm_data_r);
            end
            %title(var_list{iDOF}, 'Interpreter', 'none')
%             if any(iDOF == FirstCol)
%                 ylabel([y_label ' of the right leg'])
%             end
        end
        m = m+1;
        %add line indicating right toe off (average of all trials)
        x_val = ev.cycle.right.footOff-ev.cycle.right.start;
        x_val_norm_r = x_val_norm_r + (x_val / (ev.cycle.right.end - ev.cycle.right.start)*100); % time normalize
    end
    % plot left leg
    if isfield(ev.cycle,'left')
        for iDOF = length(var_list)/2+1:length(var_list)
            %             axes(ax(iDOF)); hold on; grid on
%             switch var_list{iDOF}
%                 case 'glut_max'
%                     glut_max = 1;
%             end
            raw_data_l = data.(var_list{iDOF})(L_GC{1}(1):L_GC{1}(2));
            time_norm_data_l = TimeNorm(raw_data_l,fs);
            if length(ev.cycle.left.start) > 1
                for k = 2:length(ev.cycle.left.start)
                    raw_data_l = TimeNorm(data.(var_list{iDOF})(L_GC{k}(1):L_GC{k}(2)),fs);
                    time_norm_data_l = time_norm_data_l + raw_data_l;
                end
                time_norm_data_l = time_norm_data_l./length(ev.cycle.left.start);
            end
            if BW_norm == 1
                weight_norm_data_l = time_norm_data_l./(subj_mass*9.81);
                data_l{iDOF-length(var_list)/2}(n,:) = weight_norm_data_l;
                %plot(weight_norm_data_l);
            else
                data_l{iDOF-length(var_list)/2}(n,:) = time_norm_data_l;
                %plot(time_norm_data_l);
            end
            %             raw_data = kin.(angles_var_list{iDOF})(L_GC(1):L_GC(2));
            %             time_norm_data = TimeNorm(raw_data,fs);
            %title(var_list{iDOF}, 'Interpreter', 'none')
%             if any(iDOF == FirstCol)
%                 ylabel([y_label ' of the left leg'])
%             end
        end
        n = n+1;
        % add line indicating left toe off (average of all trials)
        x_val = ev.cycle.left.footOff-ev.cycle.left.start;
        x_val_norm_l = x_val_norm_l + (x_val / (ev.cycle.left.end - ev.cycle.left.start)*100); % time normalize
    end
end

% figure('Name',f)
% if isfield(ev.cycle,'right')
% x_val_norm_r = x_val_norm_r/nb_r;
% for k = 1:length(ax)/2
%     axes(ax(k));
%     plotVert(x_val_norm_r);
% end
% end

% if isfield(ev.cycle,'left')
% x_val_norm_l = x_val_norm_l/nb_l;
% for k = length(ax)/2+1:length(ax)
%     axes(ax(k));
%     plotVert(x_val_norm_l);
% end
% end

% tight_subplot_ticks (ax,LastRow,0)
% 
% mmfn
end