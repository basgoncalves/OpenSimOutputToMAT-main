function [Markers_coordinates,Marker_labels] = conver_struct_2_double(trc_struct)
Marker_labels = fields(trc_struct);
Markers_coordinates = [];
for i = 1:length(Marker_labels)                    % convert trc struct into double (data) and cell (lables)
    field_data = trc_struct.(Marker_labels{i});
    for col = 1:size(field_data,2)
        Markers_coordinates(:,end+1) = field_data(:,col);
    end
end

