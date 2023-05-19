%% script combining conversion from c3d format to trc+mot and the TorsionTool
clc
clear all
close all

subj_mass = 32.9;

%% muscle forces
% SO = load_sto_file('C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\C3D\dynamic03\deformed_model-scaled_StaticOptimization_force.sto');
act = load_sto_file('C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\dyn_03\dyn_03\deformed_model-scaled_StaticOptimization_activation.sto');data = act;


ev = importdata('C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\C3D\dynamic03\settings.mat');
if isfield(ev.cycle,'left')
    L_GC = ev.cycle.left.start(1);
    L_GC(2) = ev.cycle.left.end(1);
end
if isfield(ev.cycle,'right')
    R_GC = ev.cycle.right.start(1);
    R_GC(2) = ev.cycle.right.end(1);
end
% str = data_SO_BW(hs1:to,:); % extract data of one stride

time = data.time; %time steps
% time_str = time(hs1:to,:);
if exist('R_GC')
    time_str = time(R_GC(1):R_GC(2),:);
    hs1 = R_GC(1);
    hs2 = R_GC(2);
    to_val = ev.cycle.right.footOff-ev.cycle.right.start; % toe off
    to_val_norm = to_val / (ev.cycle.right.end - ev.cycle.right.start)*100; % time normalize
    s_name = 'right';
    suff_var = '_r';
    suff_name = 'right';
elseif exist('L_GC')
    time_str = time(L_GC(1):L_GC(2),:);
    hs1 = L_GC(1);
    hs2 = L_GC(2);
    to_val = ev.cycle.left.footOff(1)-ev.cycle.left.start(1); % toe off
    to_val_norm = to_val/(ev.cycle.left.end(1) - ev.cycle.left.start(1))*100; % time normalize
    s_name = 'left';
    suff_var = '_l';
    suff_name = 'left';
end

fs = 1 / (time(2)-time(1));

% muscles_var_list = {{'glut_max1', 'glut_max2', 'glut_max3'}, {'med_gas'},  {'lat_gas'}, {'bifemlh'},{'rect_fem'},{'soleus'};...
%                     'gluteus maximus', 'gastrocnemius mediales',  'gastrocnemius lateralis', 'biceps femoris long head','rectus femoris','soleus'};
muscles_var_list = {{'tib_ant'}, {'med_gas'},  {'lat_gas'}, {'flex_hal'},{'rect_fem'},{'soleus'};...
                    'tibialis anterior', 'gastrocnemius mediales',  'gastrocnemius lateralis', 'flexor hallucis','rectus femoris','soleus'};
[ax, pos,FirstCol,LastRow,LastCol]  = tight_subplotBG(length(muscles_var_list),0);

%left leg
% suff_var = '_l';
% suff_name = 'left';
for i = 1: length(muscles_var_list)
    axes(ax(i)); hold on; grid on
        for j=1:length(muscles_var_list{1,i})
            data = data.([muscles_var_list{1,i}{j} suff_var])(hs1:hs2,:);
%             data = SO.([muscles_var_list{1,i}{j} suff_var])(hs1:hs2,:)/(subj_mass*9.81);
            data_time_norm = TimeNorm (data,fs);
            plot(data_time_norm)
            plotVert(to_val_norm)
        end
    title(muscles_var_list{1,i}, 'Interpreter', 'none')
    if any(i == FirstCol)
%         ylabel('muscle force [N/BW]')
        ylabel('muscle activation')
    end
end

tight_subplot_ticks(ax,LastRow,0)
mmfn



%% modify max isometric force
import org.opensim.modeling.*
ModelIn = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\Model\FINAL_PERSONALISEDTORSIONS_scaled_final.osim';
model = Model(ModelIn);
muscles = model.getMuscles();
nMuscles = muscles.getSize();

for ii = 0:nMuscles-1
    current_max = muscles.get(ii).getMaxIsometricForce;
    muscles.get(ii).setMaxIsometricForce(current_max *1.5);
end

% Write the model to a new file
ModelOut = strrep(ModelIn, '.osim','_new_isom.osim');
model.print(ModelOut)
%% velocity & joint moments
IK = importdata('C:\Users\Kira\Nextcloud\Documents\MA\Daten\1\pre\dyn_03\IK_dyn03.mot');
pelvis_vx = IK.data(:,5);
v_pelvis = calcVelocity(pelvis_vx,200);% walking velocity

ID = load_sto_file('C:\Users\Kira\Nextcloud\Documents\MA\Daten\1\pre\dyn_03\ID_dyn03.sto');
% data_ID = ID.data;
fr_start = 257; %nb of frames before measurement starts --> to be substracted
hs1 = 515; % frame where first heel strike occurs
hs1 = hs1 - fr_start;
hs2 = 671; %frame where second heel strike occurs
hs2 = hs2 - fr_start;
to = 599; % toe off
to = to - fr_start;
% str = data_ID(hs1:to,:); % extract data of one stride
% stp = data_ID(hs1:hs2,:); % extract data of one step
% stp = stp/(32.9*9-81); %normalization
time = ID.time; %time steps
time_str = time(hs1:to,:);
time_stp = time(hs1:hs2,:);
time_norm = 0:100/(length(time_stp)-1):100;

[hip_flex] = moments_plot(ID.hip_flexion_r_moment(hs1:hs2,:)/(32.9*9-81), time_norm, 'hip flexion', hs1, to, 'left');
[knee_angle] = moments_plot(ID.knee_angle_r_moment/(32.9*9-81), time_norm, 'knee angle', hs1, to, 'left');
[ankle_angle] = moments_plot(ID.ankle_angle_r_moment/(32.9*9-81), time_norm, 'ankle angle', hs1, to, 'left');

%% Hans' reference

% moments
ID_Hans = load_sto_file('C:\Users\Kira\Nextcloud\Documents\MA\Daten\reference Hans\OSoutputExample\NSAplus15AntePlus15_TD_FPA\ID\NSAplus15AntePlus15_TD_FPA_ID.sto');
time_stp = ID_Hans.time; %time steps
time_norm = 0:100/(length(time_stp)-1):100; % normalize time

[hip_flex] = moments_plot(ID_Hans.hip_flexion_r_moment/(32.9*9-81), time_norm, 'reference: hip flexion', 1, []);
[knee_angle] = moments_plot(ID_Hans.knee_angle_r_moment/(32.9*9-81), time_norm, 'reference: knee angle', 1, []);
[ankle_angle] = moments_plot(ID_Hans.ankle_angle_r_moment/(32.9*9-81), time_norm, 'reference: ankle angle', 1, []);

% muscle forces
SO_Hans = load_sto_file('C:\Users\Kira\Nextcloud\Documents\MA\Daten\reference Hans\OSoutputExample\NSAplus15AntePlus15_TD_FPA\SO\NSAplus15AntePlus15_StaticOptimization_force.sto');
% data_SO_BW = data_SO/(32.9*9.81); % muscle forces per body weight
% stp = data_SO_BW(hs1:hs2,:); % extract data of one step
% time_stp = SO_Hans.time;
% time_norm = 0:100/(length(time_stp)-1):100;

[glut_max_stp] = muscle_force_plot([SO_Hans.glut_max1_r/(32.9*9.81), SO_Hans.glut_max2_r/(32.9*9.81), SO_Hans.glut_max3_r/(32.9*9.81)], time_norm, 1:3, 'gluteus maximus (reference)', [], []);
[gastr_med_stp] = muscle_force_plot(SO_Hans.med_gas_r/(32.9*9.81), time_norm, 1, 'gastrocnemius medialis (reference)', [], []);
[gastr_lat_stp] = muscle_force_plot(SO_Hans.lat_gas_r/(32.9*9.81), time_norm, 1, 'gastrocnemius lateralis (reference)', [], []);
[biceps_long_stp] = muscle_force_plot(SO_Hans.bifemlh_r/(32.9*9.81), time_norm, 1, 'biceps femoris long head (reference)', [], []);
[rect_fem_stp] = muscle_force_plot(SO_Hans.rect_fem_r/(32.9*9.81), time_norm, 1, 'rectus femoris (reference)', [], []);
[sol_stp] = muscle_force_plot(SO_Hans.soleus_r/(32.9*9.81), time_norm, 1, 'soleus (reference)', [], []);

%% compare with and without grf2
clc
clear all
close all

%with grf2
hs1 = 245;
to = 450;
hs2 = 512;
fr_start = 121; %nb of frames before measurement starts --> to be substracted
hs1 = hs1 - fr_start;
hs2 = hs2 - fr_start;
to = to - fr_start;
ID2 = load_sto_file('C:\Users\Kira\Nextcloud\Documents\MA\Daten\1\pre\dyn_03\ID_dyn03.sto'); % with grf2
ID0 = load_sto_file('C:\Users\Kira\Nextcloud\Documents\MA\Daten\1\pre\dyn_03\ID_no_grf2.sto'); % without grf2
time = ID2.time; %time steps
time_stp = time(hs1:hs2,:);
time_norm = 0:100/(length(time_stp)-1):100;

figure(1)
subplot(3,1,1)
plot(ID2.time, ID2.pelvis_tilt_moment)
hold on
plot(ID2.time, ID0.pelvis_tilt_moment)
title('pelvis tilt moment with grf2')

subplot(3,1,2)
plot(ID2.time, ID2.pelvis_list_moment)
hold on
plot(ID2.time, ID0.pelvis_list_moment)
title('pelvis list moment with grf2')

subplot(3,1,3)
plot(ID2.time, ID2.pelvis_rotation_moment)
hold on
title('pelvis rotation moment with grf2')
plot(ID2.time, ID0.pelvis_rotation_moment)
title('pelvis rotation moment without grf2')
legend('ID2', 'ID0')

figure(2)
subplot(3,1,1)
plot(ID2.time, ID2.pelvis_tx_force)
hold on
plot(ID2.time, ID0.pelvis_tx_force)
title('pelvis tilt moment with grf2')
legend('ID0', 'ID2')

subplot(3,1,2)
plot(ID2.time, ID2.pelvis_ty_force)
hold on
plot(ID2.time, ID0.pelvis_ty_force)
title('pelvis list moment with grf2')

subplot(3,1,3)
plot(ID2.time, ID2.pelvis_tz_force)
hold on
plot(ID2.time, ID0.pelvis_tz_force)
title('pelvis rotation moment with grf2')

%% delete 
data = btk_loadc3d('C:\Users\Kira\Nextcloud\Documents\MA\Daten\1\pre\dyn_03\dynamic03.c3d');
% data = btkRemovePoint(data, data.Markers.