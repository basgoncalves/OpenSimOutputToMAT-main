function  plot_data_stride(data_path, subj_mass, ev_path, var_list, norm_BW, label_y)
% plot IK/SO results of one stride
data = load_sto_file(data_path);

ev = importdata(ev_path);
if isfield(ev.cycle,'left')
    L_GC = ev.cycle.left.start(1);
    L_GC(2) = ev.cycle.left.end(1);
end
if isfield(ev.cycle,'right')
    R_GC = ev.cycle.right.start(1);
    R_GC(2) = ev.cycle.right.end(1);
end

time = data.time; %time steps
if exist('R_GC')
%     time_str = time(R_GC(1):R_GC(2),:);
    hs1 = R_GC(1);
    hs2 = R_GC(2);
    to_val = ev.cycle.right.footOff-ev.cycle.right.start; % toe off
    to_val_norm = to_val / (ev.cycle.right.end - ev.cycle.right.start)*100; % time normalize
    s_name = 'right';
    suff_var = '_r';
elseif exist('L_GC')
%     time_str = time(L_GC(1):L_GC(2),:);
    hs1 = L_GC(1);
    hs2 = L_GC(2);
    to_val = ev.cycle.left.footOff(1)-ev.cycle.left.start(1); % toe off
    to_val_norm = to_val/(ev.cycle.left.end(1) - ev.cycle.left.start(1))*100; % time normalize
    s_name = 'left';
    suff_var = '_l';
end

fs = 1 / (time(2)-time(1));

[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(length(var_list),0);

for i = 1: length(var_list)
    axes(ax(i)); hold on; grid on
        for j=1:length(var_list{1,i})
            if norm_BW == 1
                data_muscle = data.([var_list{1,i}(j) suff_var])(hs1:hs2,:)/(subj_mass*9.81); % for SO
            else
                data_muscle = data.([var_list{1,i} suff_var])(hs1:hs2,:);
            end
            
            data_time_norm = TimeNorm (data_muscle, fs);
            plot(data_time_norm)
            plotVert(to_val_norm)
        end
    title(var_list{1,i}, 'Interpreter', 'none')
    if any(i == FirstCol)
        ylabel([label_y s_name ' leg'])
    end
end

tight_subplot_ticks(ax,LastRow,0)
mmfn
end