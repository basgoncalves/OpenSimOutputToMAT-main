function [dist_max, max_err, loc_err] = IK_max_marker_error(subj,id, surg, markers_torso_list, markers_legs_list)
markers_e_filepath = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\C3D\' id '\marker_experimental.trc'];
markers_v_filepath = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\C3D\' id '\_ik_model_marker_locations.sto'];
modelpath = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_PERSONALISEDTORSIONS_scaled_final.osim'];
setupxmlpath = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\IK_setup_' subj '.xml'];
save_folder = fileparts(markers_e_filepath);
outputPath = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\C3D\' id '\'];
markers_e = load_marker_trc(markers_e_filepath);
markers_list = [markers_torso_list, markers_legs_list];

for j = 1:length(markers_list)
    name_e_x = [markers_list{j}, '_X'];
    name_e_y = [markers_list{j}, '_Y'];
    name_e_z = [markers_list{j}, '_Z'];
    for k = 1:length(markers_e.(name_e_x))
        switch markers_e.(name_e_x){k}
            case 'nan'
                markers_e.(name_e_x){k} =markers_e.(name_e_x){k-1};
        end
        switch markers_e.(name_e_y){k}
            case 'nan'
                markers_e.(name_e_y){k} =markers_e.(name_e_y){k-1};
        end
        switch markers_e.(name_e_z){k}
            case 'nan'
                markers_e.(name_e_z){k} =markers_e.(name_e_z){k-1};
        end
    end
end

import org.opensim.modeling.*;
ikTool = InverseKinematicsTool(setupxmlpath);
osimModel = Model(modelpath);
ikTool.setModel(osimModel);
ikTool.setMarkerDataFileName(markers_e_filepath);
ikTool.setStartTime(markers_e.Time{1});
ikTool.setEndTime(markers_e.Time{end});
ikTool.setOutputMotionFileName(fullfile(outputPath, 'IK_Results.mot'));
ikTool.set_results_directory(outputPath)
ikTool.print(fullfile(outputPath, 'ikSettings.xml'));
ikTool.set_report_marker_locations(true);
ikTool.run();

% calculate marker errors
markers_v = load_sto_file(markers_v_filepath);
% markers_list = {'T10','C7','LASI', 'RASI', 'SACR', 'LT1', 'LT2', 'LT3', 'LS1', 'LS2', 'LS3', 'LHEE', 'LTOE', 'RT1', 'RT2', 'RT3', 'RS1', 'RS2', 'RS3', 'RHEE', 'RTOE'};
% legs_markers_list = markers_list(4:end);
max_err = 0; %magnitude of maximum marker error
loc_err = '';

for i=1:length(markers_list)
    name_e_x = [markers_list{i}, '_X'];
    name_e_y = [markers_list{i}, '_Y'];
    name_e_z = [markers_list{i}, '_Z'];
    name_v_x = [markers_list{i}, '_tx'];
    name_v_y = [markers_list{i}, '_ty'];
    name_v_z = [markers_list{i}, '_tz'];
    m_e.(markers_list{i}) = cell2mat([markers_e.(name_e_x), markers_e.(name_e_y), markers_e.(name_e_z)])./1000;
    m_v.(markers_list{i}) = [markers_v.(name_v_x), markers_v.(name_v_y), markers_v.(name_v_z)];
    dist.(markers_list{i}) = dist_markers(m_e.(markers_list{i}), m_v.(markers_list{i}));
    dist_max.(markers_list{i}) = max(dist.(markers_list{i}));

    if ismember(markers_list(i),markers_legs_list) && max(dist.(markers_list{i})) > max_err
        max_err = max(dist.(markers_list{i}));
        loc_err = markers_list{i};
    end
end
end