% rename virtual markers
% activate_msk_modelling

% load marker file
path = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD07\Model\FINAL_MARKERSET.xml';
v_markers = xml_read(path);

% rename markers
for i = 1 : length(v_markers.MarkerSet.objects.Marker)
    switch v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name
        case 'RT1'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'RTHI1';
        case 'RT2'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'RTHI2';
        case 'RT3'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'RTHI3';
        case 'RS1'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'RTIB1';
        case 'RS2'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'RTIB2';
        case 'RS3'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'RTIB3';
        case 'LT1'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'LTHI1';
        case 'LT2'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'LTHI2';
        case 'LT3'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'LTHI3';
        case 'LS1'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'LTIB1';
        case 'LS2'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'LTIB2';
        case 'LS3'
            v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name = 'LTIB3';
       
    end
end

% save new marker set
path_output = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD07_new_virtual_markers.xml';
root = 'OpenSimDocument';                                                        
Pref = struct;
Pref.StructItem = false;
Pref.CellItem = false;
v_markers = ConvertLogicToString(v_markers);
path_file_output = path_output;
xml_write(path_file_output, v_markers, root,Pref);