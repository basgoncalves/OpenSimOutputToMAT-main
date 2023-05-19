close all; clear all; clc;
fname_rTib = 'C:\Users\Hans\Documents\UNI_WIEN\KranzlCooperation\HansGriffithData_TD\MRI\slicer\TD01/RightTib_TD01.mrk.json';
fname_rFem = 'C:\Users\Hans\Documents\UNI_WIEN\KranzlCooperation\HansGriffithData_TD\MRI\slicer\TD01/RightFem_TD01.mrk.json';
fname_lTib = 'C:\Users\Hans\Documents\UNI_WIEN\KranzlCooperation\HansGriffithData_TD\MRI\slicer\TD01/LeftTib_TD01.mrk.json';
fname_lFem = 'C:\Users\Hans\Documents\UNI_WIEN\KranzlCooperation\HansGriffithData_TD\MRI\slicer\TD01/LeftFem_TD01.mrk.json';

%% get data

fid = fopen(fname_rTib);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
data_rTib = jsondecode(str);
for i = 1 : length(data_rTib.markups.controlPoints)
    switch data_rTib.markups.controlPoints(i).label
        case 'MT'
            mtI = i;
        case 'LT'
            ltI = i;
        case 'MMAL'
            mmalI = i;
        case 'LMAL'
            lmalI = i;        
    end
end

fid = fopen(fname_lTib);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
data_lTib = jsondecode(str);
for i = 1 : length(data_lTib.markups.controlPoints)
    switch data_lTib.markups.controlPoints(i).label
        case 'MT'
            mtI = i;
        case 'LT'
            ltI = i;
        case 'MMAL'
            mmalI = i;
        case 'LMAL'
            lmalI = i;        
    end
end

fid = fopen(fname_rFem);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
data_rFem = jsondecode(str);
for i = 1 : length(data_rFem.markups.controlPoints)
    switch data_rFem.markups.controlPoints(i).label
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

fid = fopen(fname_lFem);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
data_lFem = jsondecode(str);
for i = 1 : length(data_lFem.markups.controlPoints)
    switch data_lFem.markups.controlPoints(i).label
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
%% Rajagopal base model measurements
Rajagopal.Markers.L_HJC = [-0.056276 0.85151 -0.07726];
Rajagopal.Markers.R_HJC = [-0.056276 0.85151 0.07726];
Rajagopal.Markers.R_GT = [-0.0788038 0.80126 0.13722];
Rajagopal.Markers.R_PS = [-0.0662013 0.780346 0.116023];
Rajagopal.Markers.R_DS = [-0.0467056 0.439209 0.07726];
Rajagopal.Markers.R_MC = [-0.0790706 0.439788 0.0532338];
Rajagopal.Markers.R_LC = [-0.0790706 0.439788 0.100893];

Rajagopal.Markers.R_MT = [-0.0825586 0.408608 0.0648141];
Rajagopal.Markers.R_LT = [-0.0806764 0.408608 0.0978893];
Rajagopal.Markers.R_MMAL = [-0.0557446 0.0660587 0.048199];
Rajagopal.Markers.R_LMAL = [-0.0710674 0.0613259 0.113142];

dist_curr = Rajagopal.Markers.L_HJC - Rajagopal.Markers.R_HJC;
Rajagopal.interHJC = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
clear dist_curr;
dist_curr = Rajagopal.Markers.R_HJC - Rajagopal.Markers.R_LC;
Rajagopal.femLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
clear dist_curr;
dist_curr = Rajagopal.Markers.R_LT - Rajagopal.Markers.R_LMAL;
Rajagopal.tibLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
clear dist_curr;
%% calculate scale factors

dist_curr = ((data_lFem.markups.controlPoints(hjcI).position - data_rFem.markups.controlPoints(hjcI).position)/1000)';
Subj.interHJC = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scalePelvisWidth = Subj.interHJC / Rajagopal.interHJC ;
clear dist_curr;
dist_curr = ((data_rFem.markups.controlPoints(hjcI).position - data_rFem.markups.controlPoints(lcI).position)/1000)';
Subj.femLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_R_FemLength = Subj.femLength / Rajagopal.femLength ;
clear dist_curr;
dist_curr = ((data_rTib.markups.controlPoints(ltI).position - data_rTib.markups.controlPoints(lmalI).position)/1000)';
Subj.tibLength = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_R_TibLength = Subj.tibLength / Rajagopal.tibLength ;
clear dist_curr;
% left leg
dist_curr = ((data_lFem.markups.controlPoints(hjcI).position - data_lFem.markups.controlPoints(lcI).position)/1000)';
Subj.femLength_L = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_L_FemLength = Subj.femLength_L / Rajagopal.femLength ;
clear dist_curr;
dist_curr = ((data_lTib.markups.controlPoints(ltI).position - data_lTib.markups.controlPoints(lmalI).position)/1000)';
Subj.tibLength_L = sqrt(dist_curr(1,1)^2 + dist_curr(1,2)^2 + dist_curr(1,3)^2);
Subj.scale_L_TibLength = Subj.tibLength_L / Rajagopal.tibLength ;
clear dist_curr;


