function [scalePelvis, scaleTibR, scaleTibL, scaleFemR, scaleFemL] = scale_template(subj,subj_name, subj_mass, surg)

switch surg
    case '\post'
        v_marker_set = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_MARKERSET.xml']; % virtual marker set after Torsion Tool
    case ''
        v_marker_set = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_MARKERSET_MRI.xml']; % virtual marker set after Torsion Tool
    case '\pre'
        v_marker_set = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_MARKERSET_MRI.xml']; % virtual marker set after Torsion Tool
end
% v_marker_set = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_MARKERSET_MRI.xml']; % virtual marker set after Torsion Tool
fname_R = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\Markups\' subj '_R_Fiducials.mrk.json'];
fname_L = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\Markups\' subj '_L_Fiducials.mrk.json'];
path_output = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Scale_template_' subj '.xml'];
path_markerset = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_MARKERSET.xml'];
output_motion_path = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_PERSONALISEDTORSIONS_scaled_final.mot'];
% output_model_path

%% modify scaling template
sc_template = xml_read('C:\Users\Balu\Nextcloud\Documents\MA\Daten\ScaleTool_template.xml');
sc_template.ScaleTool.mass = subj_mass;
sc_template.ScaleTool.GenericModelMaker.marker_set_file = path_markerset;
sc_template.ScaleTool.ATTRIBUTE.name = subj_name;
sc_template.ScaleTool.MarkerPlacer.output_motion_file = output_motion_path;

v_markers = xml_read(v_marker_set);

switch surg
    case ''
    % calculate scale factors
    Subj = CalculateScaleFactors_MRIKira(fname_R,fname_L, v_markers);
    scalePelvis = Subj.scalePelvisWidth;
    scaleTibR = Subj.scale_R_TibLength;
    scaleTibL = Subj.scale_L_TibLength;
    scaleFemR = Subj.scale_R_FemLength;
    scaleFemL = Subj.scale_L_FemLength;

    % change scale factors in template
    for i = 1 : length(sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale)
        switch sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).segment
            case 'pelvis'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scalePelvis, scalePelvis, scalePelvis];
            case 'femur_r'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleFemR, scaleFemR, scaleFemR];
            case 'femur_l'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleFemL, scaleFemL, scaleFemL];
            case 'tibia_r'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleTibR, scaleTibR, scaleTibR];
            case 'tibia_l'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleTibL, scaleTibL, scaleTibL];
        end
    end
    case '\pre'
    % calculate scale factors
    Subj = CalculateScaleFactors_MRIKira(fname_R,fname_L, v_markers);
    scalePelvis = Subj.scalePelvisWidth;
    scaleTibR = Subj.scale_R_TibLength;
    scaleTibL = Subj.scale_L_TibLength;
    scaleFemR = Subj.scale_R_FemLength;
    scaleFemL = Subj.scale_L_FemLength;

    % change scale factors in template
    for i = 1 : length(sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale)
        switch sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).segment
            case 'pelvis'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scalePelvis, scalePelvis, scalePelvis];
            case 'femur_r'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleFemR, scaleFemR, scaleFemR];
            case 'femur_l'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleFemL, scaleFemL, scaleFemL];
            case 'tibia_r'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleTibR, scaleTibR, scaleTibR];
            case 'tibia_l'
                sc_template.ScaleTool.ModelScaler.ScaleSet.objects.Scale(i).scales = [scaleTibL, scaleTibL, scaleTibL];
        end
    end
end



%% save new scaling template
root = 'OpenSimDocument';                                                        
Pref = struct;
Pref.StructItem = false;
Pref.CellItem = false;
sc_template = ConvertLogicToString(sc_template);
path_file_output = path_output;
xml_write(path_file_output, sc_template, root,Pref);
end