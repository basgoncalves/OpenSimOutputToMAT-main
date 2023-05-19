function scale_template_no_mri(subj,subj_name, subj_mass, surg)


path_output = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Scale_template_' subj '.xml'];
path_markerset = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_MARKERSET.xml'];
output_motion_path = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_PERSONALISEDTORSIONS_motion.mot'];
output_model_path = ['C:\Users\Balu\Nextcloud\Documents\MA\Daten\' subj surg '\Model\FINAL_PERSONALISEDTORSIONS_scaled.osim'];

%% modify scaling template
sc_template = xml_read('C:\Users\Balu\Nextcloud\Documents\MA\Daten\ScaleTool_template_no_mri.xml');
sc_template.ScaleTool.mass = subj_mass;
sc_template.ScaleTool.GenericModelMaker.marker_set_file = path_markerset;
sc_template.ScaleTool.ATTRIBUTE.name = subj_name;
sc_template.ScaleTool.MarkerPlacer.output_motion_file = output_motion_path;
sc_template.ScaleTool.MarkerPlacer.output_model_file = output_model_path;

%% save new scaling template
root = 'OpenSimDocument';                                                        
Pref = struct;
Pref.StructItem = false;
Pref.CellItem = false;
sc_template = ConvertLogicToString(sc_template);
path_file_output = path_output;
xml_write(path_file_output, sc_template, root,Pref);
end