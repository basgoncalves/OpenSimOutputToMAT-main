function [Subj] = CalculateScaleFactors_MRIKira(fname_R,fname_L, v_markers)

%% get data
% switch surg
%     case ''

[data_L,ltI_l,lmalI_l,hjcI_l, lcI_l] = find_i_markups(fname_L);
[data_R,ltI,lmalI,hjcI, lcI] = find_i_markups(fname_R);
% data_R.markups.controlPoints(hjcI).position
% data_R.markups.controlPoints(lcI).position
% data_L.markups.controlPoints(hjcI_l).position
% data_L.markups.controlPoints(lcI_l).position
%% OS2392 base model measurements after Torsion Tool
% data.markups.controlPoints(hjcI).position = [-0.070699999999999999; 0.88389999999999991; 0.083500000000000005];

for i = 1 : length(v_markers.MarkerSet.objects.Marker)
    switch v_markers.MarkerSet.objects.Marker(i).ATTRIBUTE.name  
%         case 'RHJC'
%             OS2392.RHJC = v_markers.MarkerSet.objects.Marker(i).location;
%         case 'LHJC'
%             OS2392.LHJC = v_markers.MarkerSet.objects.Marker(i).location;
        case 'LLT'
            OS2392.Markers.L_LT = v_markers.MarkerSet.objects.Marker(i).location;
        case 'LLC'
            OS2392.Markers.L_LC = v_markers.MarkerSet.objects.Marker(i).location;
        case 'LANK'
            OS2392.Markers.L_LMAL = v_markers.MarkerSet.objects.Marker(i).location;
        case 'LHJC_g'
            OS2392.Markers.L_HJC_g = v_markers.MarkerSet.objects.Marker(i).location;
        case 'LHJC'
            OS2392.Markers.L_HJC = v_markers.MarkerSet.objects.Marker(i).location;
        case 'RLT'
            OS2392.Markers.R_LT = v_markers.MarkerSet.objects.Marker(i).location;
        case 'RLC'
            OS2392.Markers.R_LC = v_markers.MarkerSet.objects.Marker(i).location;
        case 'RANK'
            OS2392.Markers.R_LMAL = v_markers.MarkerSet.objects.Marker(i).location;
        case 'RHJC_g'
            OS2392.Markers.R_HJC_g = v_markers.MarkerSet.objects.Marker(i).location;
        case 'RHJC'
            OS2392.Markers.R_HJC = v_markers.MarkerSet.objects.Marker(i).location;
    end
end

if isfield(OS2392.Markers, 'R_LC') == 0
    OS2392.Markers.R_LT = OS2392.Markers.L_LT; OS2392.Markers.R_LT(3) = -OS2392.Markers.R_LT(3);
    OS2392.Markers.R_LC = OS2392.Markers.L_LC; OS2392.Markers.R_LC(3) = -OS2392.Markers.R_LC(3);
end

dist_curr = OS2392.Markers.L_HJC_g - OS2392.Markers.R_HJC_g;
OS2392.interHJC = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2); % distance between HJC = width of pelvis
clear dist_curr;
dist_curr = OS2392.Markers.R_HJC - OS2392.Markers.R_LC;
OS2392.R_femLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2); % length of right femur
clear dist_curr;
dist_curr = OS2392.Markers.R_LT - OS2392.Markers.R_LMAL;
OS2392.R_tibLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2); % length of left tibia
clear dist_curr;
dist_curr = OS2392.Markers.L_HJC - OS2392.Markers.L_LC;
OS2392.L_femLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);% length of left femur
clear dist_curr;
dist_curr = OS2392.Markers.L_LT - OS2392.Markers.L_LMAL;
OS2392.L_tibLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2); % length of left tibia
clear dist_curr;
%% calculate scale factors

dist_curr = ((data_L.markups.controlPoints(hjcI_l).position - data_R.markups.controlPoints(hjcI).position)/1000)';
Subj.interHJC = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scalePelvisWidth = Subj.interHJC / OS2392.interHJC;
clear dist_curr;
% right leg
dist_curr = ((data_R.markups.controlPoints(hjcI).position - data_R.markups.controlPoints(lcI).position)/1000)';
Subj.femLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_R_FemLength = Subj.femLength / OS2392.R_femLength ;
clear dist_curr;
dist_curr = ((data_R.markups.controlPoints(ltI).position - data_R.markups.controlPoints(lmalI).position)/1000)';
Subj.tibLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_R_TibLength = Subj.tibLength / OS2392.R_tibLength ;
clear dist_curr;
% left leg
dist_curr = ((data_L.markups.controlPoints(hjcI_l).position - data_L.markups.controlPoints(lcI_l).position)/1000)';
Subj.femLength_L = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_L_FemLength = Subj.femLength_L / OS2392.L_femLength;
clear dist_curr;
dist_curr = ((data_L.markups.controlPoints(ltI_l).position - data_L.markups.controlPoints(lmalI_l).position)/1000)';
Subj.tibLength_L = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_L_TibLength = Subj.tibLength_L / OS2392.L_tibLength
clear dist_curr;


function [data,ltI,lmalI,hjcI, lcI] = find_i_markups(fname)

fid = fopen(fname);
raw = fread(fid,inf);   
str = char(raw');
fclose(fid);
data = jsondecode(str);

for i = 1 : length(data.markups.controlPoints)
    switch data.markups.controlPoints(i).label
        case 'MT'
            mtI = i;
        case 'LT'
            ltI = i;
        case 'MMAL'
            mmalI = i;
        case 'LMAL'
            lmalI = i;
        case 'HJC'
            hjcI = i;
        case 'GT'
            gtI = i;
        case 'PS'
            psI = i;
        case 'DS'
            dsI = i;
        case 'LC'
            lcI = i;
        case 'MC'
            mcI = i;
    end
end
