function replace_marker_names_Hans(trc_path,new_trc_path)

trialpath_trc_template = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD06\C3D\TD06 Cal 01\marker_experimental_with_sacrum.trc';

trc_template = load_trc_file(trialpath_trc_template);
trc_to_edit = load_trc_file(trc_path);


marker_names = {'MKNE','TH1', 'TH2','TH3','MMA','TB1','TB2','TB3', 'HEE'};
replace_names = {'KNEM','THI1','THI2','THI3','ANKM','TIB1','TIB2','TIB3', 'HEEL'};

marker_names = [cellfun(@(x) ['L' x], marker_names, 'UniformOutput', false) cellfun(@(x) ['R' x], marker_names, 'UniformOutput', false)];
replace_names = [cellfun(@(x) ['L' x], replace_names, 'UniformOutput', false) cellfun(@(x) ['R' x], replace_names, 'UniformOutput', false)];


for i = 1:length(marker_names)
    try
        trc_to_edit.(replace_names{i}) = trc_to_edit.(marker_names{i});
        trc_to_edit = rmfield(trc_to_edit,marker_names{i});
    catch
    end
end

try trc_to_edit.SACR = (trc_to_edit.LPSI + trc_to_edit.RPSI)./2; catch; end

try trc_to_edit.LKJC = (trc_to_edit.LKNM + trc_to_edit.LKNE)./2; catch; end
try trc_to_edit.LAJC = (trc_to_edit.LANK + trc_to_edit.LANM)./2; catch; end

try trc_to_edit.RKJC = (trc_to_edit.RKNM + trc_to_edit.RKNE)./2; catch; end
try trc_to_edit.RAJC = (trc_to_edit.RANK + trc_to_edit.RANM)./2; catch; end

flds_edit = fields(trc_to_edit);
flds_template = fields(trc_template);
not_in_edit_after = setdiff(flds_template,flds_edit);

if ~isempty(not_in_edit_after)
    warning('some markers still missing')
    disp(not_in_edit_after)
end

[Markers_coordinates,Marker_labels] = conver_struct_2_double(trc_to_edit);
Rate = 1/(trc_to_edit.Time(2) - trc_to_edit.Time(1));
if nargin < 3
    new_trc_path = strrep(trc_path,'.trc','_new.trc');
end
writetrc_os4(Markers_coordinates,Marker_labels(2:end),Rate,new_trc_path)
