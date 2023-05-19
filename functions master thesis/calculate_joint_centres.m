function trc = calculate_joint_centres(statictrcpath)
%Hip joint center computation according to Harrington et al J.Biomech 2006

new_trc_file = strrep(statictrcpath,'.trc','_HJC.trc');

% if ~was_before_april_21(new_trc_file); return, end

trc  = load_trc_file(statictrcpath);
Rate = 1/(trc.Time(2) - trc.Time(1));

if nargin > 1                                                                                                       % if includes the scaleXML_setup, remove the markers not in the template file
    trc = remove_marker_not_in_scale_setup(trc,setupScaleXml_template);
end
[trc] = add_HJC_Harrington(trc);                                                                                    % add HJC based on Harrington equations

plot_hjc_3D(trc)    % use to confirm the position of th HJC

[Markers_coordinates,Marker_labels] = conver_struct_2_double(trc);

writetrc_os4(Markers_coordinates,Marker_labels(2:end),Rate,new_trc_file)


%============================================================================================%
function trc = remove_marker_not_in_scale_setup(trc,setupScaleXml_template)
trc_markers     = fields(trc);
Scale           = xml_read(setupScaleXml_template);
MarkerSet       = Scale.ScaleTool.MarkerPlacer.IKTaskSet.objects.IKMarkerTask;
MarkerSetNames  = {};

for i = (1:length(MarkerSet))                                                                                       % convert structure of markerset to names
    MarkerSetNames{i} = MarkerSet(i).ATTRIBUTE.name;
end
MarkerSetNames = ['Time' MarkerSetNames];

for i = flip(1:length(trc_markers))                                                                                 % delete trc markers that are not in the markerset
    iName = trc_markers{i};
    if ~contains(MarkerSetNames,iName)
        trc = rmfield(trc,trc_markers{i});
        trc_markers(i) = [];
    end
end

if contains(trc_markers,'RHJC')
    return
end

%============================================================================================%
function plot_hjc_3D (trc)
% trc should be a struct resulting from "load_trc_file(trcPath)" with
% the fields RHJC, LHJC, RASI, LASI, SACR

f = figure;
hold on
labels = {'RASI','LASI','SACR','RHJC','LHJC'};
for iLab = 1:length(labels)
    point = trc.(labels{iLab})(1,:);
    scatter3(point(1), point(2), point(3), 'filled', 'MarkerFaceColor', 'r'); % RASI in red
    text(point(1)+0.05, point(2), point(3), labels{iLab} , 'Color', 'k', 'FontSize', 10); % RASI label
end

%============================================================================================%
function [trc_markers] = add_HJC_Harrington(trc_markers)
%Hip joint center computation according to Harrington et al J.Biomech 2006
%
%PW: width of pelvis (distance among ASIS)
%PD: pelvis depth = distance between mid points joining PSIS and ASIS 
%All measures are in mm
%Harrington formula:
% x= -0.24 PD-9.9
% y= -0.30PW-10.9
% z=+0.33PW+7.3
%Developed by Zimi Sawacha <zimi.sawacha@dei.unipd.it>
%Modified by Claudio Pizzolato <claudio.pizzolato@griffithuni.edu.au>

%Renamd for convenience 

LASIS = trc_markers.LASI';   %after transposition: [3xtime]
RASIS = trc_markers.RASI';
try
    SACRUM = trc_markers.SACR';
catch
    SACRUM = (trc_markers.LPSI'+trc_markers.RPSI')./2;
    trc_markers.SACR = SACRUM';
end


for t=1:size(RASIS,2)

    %Right-handed Pelvis reference system definition    
    %Global Pelvis Center position
    OP(:,t)=(LASIS(:,t)+RASIS(:,t))/2;    
    
    PROVV(:,t)=(RASIS(:,t)-SACRUM(:,t))/norm(RASIS(:,t)-SACRUM(:,t));  
    IB(:,t)=(RASIS(:,t)-LASIS(:,t))/norm(RASIS(:,t)-LASIS(:,t));    
    
    KB(:,t)=cross(IB(:,t),PROVV(:,t));                               
    KB(:,t)=KB(:,t)/norm(KB(:,t));
    
    JB(:,t)=cross(KB(:,t),IB(:,t));                               
    JB(:,t)=JB(:,t)/norm(JB(:,t));
    
    OB(:,t)=OP(:,t);
      
    %rotation+ traslation in homogeneous coordinates (4x4)
    pelvis(:,:,t)=[IB(:,t) JB(:,t) KB(:,t) OB(:,t);
                    0 0 0 1];
    
    %Trasformation into pelvis coordinate system (CS)
    OPB(:,t)=inv(pelvis(:,:,t))*[OB(:,t);1];    
       
    PW(t)=norm(RASIS(:,t)-LASIS(:,t));
    PD(t)=norm(SACRUM(:,t)-OP(:,t));
    
    %Harrington formulae (starting from pelvis center)
    diff_ap(t)=-0.24*PD(t)-9.9;
    diff_v(t)=-0.30*PW(t)-10.9;
    diff_ml(t)=0.33*PW(t)+7.3;
    
    %vector that must be subtract to OP to obtain hjc in pelvis CS
    vett_diff_pelvis_sx(:,t)=[-diff_ml(t);diff_ap(t);diff_v(t);1];
    vett_diff_pelvis_dx(:,t)=[diff_ml(t);diff_ap(t);diff_v(t);1];    
    
    %hjc in pelvis CS (4x4)
    rhjc_pelvis(:,t)=OPB(:,t)+vett_diff_pelvis_dx(:,t);  
    lhjc_pelvis(:,t)=OPB(:,t)+vett_diff_pelvis_sx(:,t);  
    
    %Transformation Local to Global
    RHJC(:,t)=pelvis(1:3,1:3,t)*[rhjc_pelvis(1:3,t)]+OB(:,t);
    LHJC(:,t)=pelvis(1:3,1:3,t)*[lhjc_pelvis(1:3,t)]+OB(:,t);
end

trc_markers.RHJC = RHJC';
trc_markers.LHJC = LHJC';
