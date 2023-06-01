
function stats_Kira()
add_repos_to_path

[Results] = create_data_struct;
S = get_subjects;

gather_data_in_struct = false;
run_stats = true;

%% organise data in struct
if gather_data_in_struct
    Analyses = get_analyses;

    for iSess = 1:length(S.sessions)
        for iAnal = 1:length(Analyses)
            curr_analysis = Analyses{iAnal};
            for iSubj = 1:length(S.subjects)

                % load paths and whether the subject/session exist
                [Logic, Path] = get_data_paths(iSubj, iSess);
                if ~Logic; continue; end

                load(Path.matResults)
                session = ['session' num2str(iSess)];

                data_struct = data.(curr_analysis).FINAL_PERSONALISEDTORSIONS_scaled_final;
                results_struct = Results.(curr_analysis).(session);

                % add all participants to the results struct and time normalise data
                [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct,iSubj);

                Results.(curr_analysis).(session) = results_struct;

            end

            % show a progress bar
            percentage = (iSubj/length(S.subjects)*100);
            LoadingBar(percentage)
        end
    end

    % add normalised values 
    Results = normalise_muscle_forces_max_isom(Results);


    % save data
    cd(get_main_dir)
    save('results.mat', 'Results')
end
%% stats
if run_stats

    load([get_main_dir fp 'results.mat'])
    Results = calculatePeaks(Results);

    load_participant_characteristics()

    %     x = Results.JRL.session1.left.peak_HCF';
    %     y = Results.JRL.session1.right.peak_HCF';
    %     [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,1,0.05);

    plotMuscleForces(Results,'Normalised')
    
end

%% -------------------------------------------------------------------------------------------------------------- %
% ---------------------------------------------------- FUCNTIONS ------------------------------------------------ %
% --------------------------------------------------------------------------------------------------------------- %
function out  = fp
out  = filesep;

% --------------------------------------------------------------------------------------------------------------- %
function LoadingBar(percentage)
% Create or update the loading bar
persistent loadingBarHandle;
if isempty(loadingBarHandle) || ~isvalid(loadingBarHandle)
    loadingBarHandle = waitbar(percentage/100, 'Progress');
else
    waitbar(percentage/100, loadingBarHandle);
end

% Close the loading bar when the completion percentage reaches 100%
if percentage >= 100
    close(loadingBarHandle);
    clear loadingBarHandle;
end

% --------------------------------------------------------------------------------------------------------------- %
function add_repos_to_path
activeFile = [mfilename('fullpath') '.m'];

OpenSimOutputToMAT_Dir = fileparts(activeFile);
addpath(genpath(OpenSimOutputToMAT_Dir))
disp([OpenSimOutputToMAT_Dir ' added to path'])

addMSK = true;

msk_modellingDir = [fileparts(OpenSimOutputToMAT_Dir) '\MSKmodelling'];
if isinpath(msk_modellingDir) && ~addMSK
    rmpath(genpath(msk_modellingDir))
    disp([msk_modellingDir ' removed from path'])
else
    addpath(genpath(msk_modellingDir))
    disp([msk_modellingDir ' added to path'])
end

% --------------------------------------------------------------------------------------------------------------- %
function out = isinpath(directory)

% Use the which function to check if the directory is in the path
fullPath = which(directory);

% Check if the fullPath is empty to determine if the directory is in the path
if isempty(fullPath)
    out = false;
else
    out = tue;
end

% --------------------------------------------------------------------------------------------------------------- %
function data_path = get_main_dir()
activeFile = [mfilename('fullpath') '.m'];
data_path = [fileparts(fileparts(activeFile)) '\Kira_MSc_data'];

% --------------------------------------------------------------------------------------------------------------- %
function S = get_subjects()
S = struct;
S.subjects = {getfolders(get_main_dir).name};
S.sessions = {'\pre', '\post', ''};

% --------------------------------------------------------------------------------------------------------------- %
function Analyses = get_analyses()
Analyses = {'IK','ID','SO','SO_Activation','JRL'};

% --------------------------------------------------------------------------------------------------------------- %
function load_participant_characteristics()
cd(get_main_dir)
load('participants_characteristics.mat')

% --------------------------------------------------------------------------------------------------------------- %
function [Results] = create_data_struct()

Results = struct;
legs = {'right','left'};
for i = 1:3
    for l = 1:numel(legs)
        session = ['session' num2str(i)];
        Results.IK.(session).(legs{l}) = struct;
        Results.ID.(session).(legs{l}) = struct;
        Results.SO.(session).(legs{l}) = struct;
        Results.SO_Activation.(session).(legs{l}) = struct;
        Results.SO_normalised.(session).(legs{l}) = struct;
        Results.JRL.(session).(legs{l}) = struct;
        Results.Mass.(session) = [];
        Results.Subject.(session) = [];
    end
end

Results = assign_nan_to_empty_columns(Results);

% --------------------------------------------------------------------------------------------------------------- %
function [variable_struct] = find_all_variables_per_analysis_type()

S = get_subjects;
Analyses = get_analyses;
variable_struct = struct;
for iSess = 1:length(S.sessions)
    for iSubj = 1:length(S.subjects)
        [Logic, Path] = get_data_paths(iSubj, iSess);
        if ~Logic; continue; end
        load(Path.matResults)
        for iAnal = 1:length(Analyses)
            curr_analysis = Analyses{iAnal};
            data_struct = data.(curr_analysis).FINAL_PERSONALISEDTORSIONS_scaled_final;
            trial_list = fields(data_struct);
            variables = fields(data_struct.(trial_list{2}));
            for iVar = 1:length(variables)
                current_variable = variables {iVar};
                variable_struct.(curr_analysis).(current_variable) = 1;
            end
        end
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function Results = assign_nan_to_empty_columns(Results)

S = get_subjects;
Analyses = get_analyses;
legs = {'right','left'};

% loop through all trials in the data folder find unique vraible names for
% each analysis (e.g. IK = pelvis_tilt, hip_flexion... | SO = add_brev_r,
% psoas_r ...)
[variable_struct] = find_all_variables_per_analysis_type;

for iSess = 1:length(S.sessions)
    for iSubj = 1:length(S.subjects)
        [Logic, Path] = get_data_paths(iSubj, iSess);
        session = ['session' num2str(iSess)];

        for iAnal = 1:length(Analyses)
            for l = 1:numel(legs)
                curr_analysis = Analyses{iAnal};
                variables = fields(variable_struct.(curr_analysis));

                for iVar = 1:length(variables)
                    current_variable = variables {iVar};
                    
                    % if data exists assign column with 1, if not assign 0
                    if Logic == 1
                        Results.(curr_analysis).(session).(legs{l}).(current_variable)(:,iSubj) = ones(101,1);
                        Results.Subject.(session)(iSubj) = 1;
                    else
                        Results.(curr_analysis).(session).(legs{l}).(current_variable)(:,iSubj) = zeros(101,1);
                        Results.Subject.(session)(iSubj) = 0;
                    end
                    
                    % when it gets to the last subject, make all ZEROS == NaN
                    if iSubj == length(S.subjects)
                        curr_var_data = Results.(curr_analysis).(session).(legs{l}).(current_variable);
                        curr_var_data(curr_var_data==0) = nan;
                        Results.(curr_analysis).(session).(legs{l}).(current_variable) = curr_var_data;
                    end
                end
            end
        end
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function [Logic, Path] = get_data_paths(iSubj,iSess)

S = get_subjects;
session_path     = [get_main_dir fp S.subjects{iSubj} fp S.sessions{iSess}];
output_folder    = [session_path fp fp 'output automization\FINAL_PERSONALISEDTORSIONS_scaled_final'];
trialNames       = {getfolders(session_path).name};

Path.main       = get_main_dir;
Path.session    = session_path;
Path.output     = output_folder;
Path.trials     = cellfun(@(c)[session_path fp c],trialNames,'uni',false);
Path.trialNames = trialNames;
Path.matResults = [fileparts(output_folder) fp 'dataStruct_ErrorScores_no_trials_removed.mat'];
Path.osimModel  = [session_path fp 'model\FINAL_PERSONALISEDTORSIONS_scaled_final.osim'];

Logic = isfolder(Path.output);
if Logic == 0
%     disp([Path.session ' does not exist'])
end

% --------------------------------------------------------------------------------------------------------------- %
function mvic = get_max_isom_force_per_muscle(model_path)
import org.opensim.modeling.*
model    = Model(model_path);
muscles  = model.getMuscles();
nMuscles = muscles.getSize();
mvic     = struct;
for ii = 0:nMuscles-1
    muscle_name = char(muscles.get(ii));
    mvic.(muscle_name) = muscles.get(ii).getMaxIsometricForce;
end

% --------------------------------------------------------------------------------------------------------------- %
function [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct,iCol)
% data_struct = a struct similar to the output of load_sto_file

% add all participants to the results struct and time normalise data
trialNames = fields(data_struct);
legs = {'left', 'right'};
if contains(trialNames{1},'mass')
    trialNames(1) = [];
end

% remove variables that contain calc
variables  = fields(data_struct.(trialNames{1}));
variables = variables(~contains(variables,'calc'));

% create structs
raw   = struct; raw.left = struct; raw.right = struct;
tnorm = struct; tnorm.left = struct; tnorm.right = struct;

% add the varaibles to the structs (e.g. pelvis_tilt, hip_flexion...)
for iVar = variables'
    raw.right.(iVar{1}) = [];
    raw.left.(iVar{1}) = [];
    tnorm.right.(iVar{1}) = [];
    tnorm.left.(iVar{1}) = [];
end

% separate right and left trials into different structs
for iTrial = trialNames'
    for iVar = variables'

        % if variable doesn't exist in one trial, make it NaN
        try
            currentTrialData = data_struct.(iTrial{1}).(iVar{1});
        catch
            currentTrialData = nan(101,1);
        end

        % time normalise
        rows = length(currentTrialData);
        t1 = data_struct.(iTrial{1}).time(1);
        t2 = data_struct.(iTrial{1}).time(2);
        fs = 1/(t2-t1);

        if contains(iTrial{1},'_right')
            raw.right.(iVar{1})(1:rows,end+1) = currentTrialData;
            tnorm.right.(iVar{1})(1:101,end+1) = TimeNorm(currentTrialData,fs);
        else
            raw.left.(iVar{1})(1:rows,end+1) = currentTrialData;
            tnorm.left.(iVar{1})(1:101,end+1) = TimeNorm(currentTrialData,fs);
        end
    end
end

% make zeros = NaNfor iTrial = trialNames'
for iVar = variables'
    raw.right.(iVar{1})(raw.right.(iVar{1})==0) = NaN;
    raw.left.(iVar{1})(raw.left.(iVar{1})==0) = NaN;
end

% mean tnorm data
for iVar = variables'
    tnorm.right.(iVar{1}) = mean(tnorm.right.(iVar{1}),2);
    tnorm.left.(iVar{1})  = mean(tnorm.left.(iVar{1}),2);
end

% add to results struct
for iVar = variables'
    for l = 1:numel(legs)
        leg = legs{l};
        rows = [1:length(tnorm.(leg).(iVar{1}))]';
        
        % if data from "iVar" is empty, make it NaN
        if ~isnan(results_struct.(leg).(iVar{1})(1,iCol))
            results_struct.(leg).(iVar{1})(rows,iCol) = tnorm.(leg).(iVar{1});
        else
            results_struct.(leg).(iVar{1})(rows,iCol) = NaN;
            disp([leg ' ' iVar{1} ' did not exist']);
        end
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function Results = normalise_muscle_forces_max_isom(Results)

S = get_subjects;
for iSess = 1:length(S.sessions)
    session = ['session' num2str(iSess)];

    % get the index of subjects that are presented in each session (i.e. that had results)
    subject_idx = Results.Subject.(session);

    % temporary variables: right and left struct to store the normalised data
    right_norm_temp = struct;
    left_norm_temp = struct;

    % number of subject indices = number of columns in each muscle per session
    for iCol = 1:length(subject_idx)
        if subject_idx(iCol) == 0
            continue
        end
        iSubj = iCol;
        [Logic, Path] = get_data_paths(iSubj, iSess);
        % get max isometric force for each muscle in the model path
        mvic = get_max_isom_force_per_muscle(Path.osimModel);

        % normalise muscle forces per muscle based on MVIC of the model
        muscle_names = fields(mvic);
        for iMusc = 1:length(muscle_names)
            current_muscle = muscle_names{iMusc};

            % try to add
            try
                right_norm_temp.(current_muscle)(:,iCol) = Results.SO.(session).right.(current_muscle)(:,iCol) / mvic.(current_muscle) * 100;
                left_norm_temp.(current_muscle)(:,iCol) = Results.SO.(session).left.(current_muscle)(:,iCol) / mvic.(current_muscle) * 100;
            catch
                right_norm_temp.(current_muscle) = [];
                left_norm_temp.(current_muscle) = [];
                right_norm_temp.(current_muscle)(:,iCol) = Results.SO.(session).right.(current_muscle)(:,iCol) / mvic.(current_muscle) * 100;
                left_norm_temp.(current_muscle)(:,iCol) = Results.SO.(session).left.(current_muscle)(:,iCol) / mvic.(current_muscle) * 100;
            end

        end
    end

    % append data to the main Results struct
    Results.SO_normalised.(session).right = right_norm_temp;
    Results.SO_normalised.(session).left = left_norm_temp;
end

% --------------------------------------------------------------------------------------------------------------- %
function Results = calculatePeaks(Results)

Legs = {'right','left'};

for i = 1:3
    for l = 1:length(Legs)
        session = ['session' num2str(i)];
        leg = Legs{l};

        % find peak HCF
        Results.JRL.(session).(leg).HCF = sum3D(Results.JRL.(session).(leg), ['hip_' leg(1) '_on_pelvis_in_pelvis']);
        Results.JRL.(session).(leg).peak_HCF = calcMaxMovAverage(Results.JRL.(session).(leg).HCF);

        % find peak KCF
        Results.JRL.(session).(leg).KCF = sum3D(Results.JRL.(session).(leg), ['knee_' leg(1) '_on_tibia_' leg(1) '_in_tibia_' leg(1)]);
        Results.JRL.(session).(leg).peak_KCF = calcMaxMovAverage(Results.JRL.(session).(leg).KCF);

        Muscles = fields(Results.SO.(session).(leg));
        for m = 1:length(Muscles)
            muscle = Muscles{m};
            Results.SO.(session).(leg).(['peak_' muscle]) = calcMaxMovAverage(Results.SO.(session).(leg).(muscle));
        end

    end
end

% --------------------------------------------------------------------------------------------------------------- %
function out = sum3D(Results_substruct, variable)
x = Results_substruct.([variable '_fx']);
y = Results_substruct.([variable '_fy']);
z = Results_substruct.([variable '_fz']);

out = sqrt(x.^2 + y.^2 + z.^2);

% --------------------------------------------------------------------------------------------------------------- %
function out = calcMaxMovAverage(time_curve)

out = max(movmean(time_curve,5));

% --------------------------------------------------------------------------------------------------------------- %
function folders = getfolders(directory, target_string, IgnoreCase)

folders = dir (directory);

% delete "../" and "./"
try
    folders(1:2) =[];
catch
    return
end

% select only elemnets that contain folders
folders = folders([folders.isdir]);

% select folders containing "target_tring"
if nargin > 1
    if IgnoreCase == 0
        folders = folders(contains({folders.name},target_string));
    else
        folders = folders(contains({folders.name},target_string,"IgnoreCase",true));
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function TimeNormalizedData = TimeNorm(Data, fs)
% TimeNorm: Perform time normalization on the input data
%   Data: Input data matrix
%   fs: Sampling frequency

TimeNormalizedData = [];

% Loop through each column of the input data
for col = 1:size(Data, 2)

    % Get the current column of data
    currentData = Data(:, col);

    % Remove any NaN values from the current column
    currentData(isnan(currentData)) = [];

    % If length of the current column is less than 3, assign NaN values to the corresponding column in the output matrix
    if length(currentData) < 3
        TimeNormalizedData(1:101, col) = NaN;
        continue; % Skip to the next iteration of the loop
    end

    % Define time values for the current trial
    timeTrial = 0 : 1/fs : size(currentData, 1)/fs;
    timeTrial(end) = []; % Remove the last time point

    % Generate normalized time values
    Tnorm = timeTrial(end)/101 : timeTrial(end)/101 : timeTrial(end);

    % Perform interpolation to obtain the time-normalized data
    TimeNormalizedData(1:101, col) = interp1(timeTrial, currentData, Tnorm)';
end






% --------------------------------------------------------------------------------------------------------------- %
% ---------------------------------------- PLOTTING FUNCTIONS --------------------------------------------------- %
% --------------------------------------------------------------------------------------------------------------- %
function [cMat,LineStyles,Marker] = getColor (pallet_name,nColors)

if nColors>163
    error('colorBG function only works for n2 =<63')
end

if nargin < 1
    cMat = convertRGB([176, 104, 16; ...
        16, 157, 176; ...
        136, 16, 176;176,...
        16, 109;31, 28, 28]);  % color scheme 2 (Bas)
else

    switch pallet_name

        case 'prism';   cMat = prism;
        case 'parula';  cMat = parula;
        case 'flag';    cMat = flag;
        case 'hsv';     cMat = hsv;
        case 'hot';     cMat = hot;
        case 'cool';    cMat = cool;
        case 'spring';  cMat = spring;
        case 'summer';  cMat = summer;
        case 'autumn';  cMat = autumn;
        case 'winter';  cMat = winter;
        case 'gray';    cMat = gray;
        case 'bone';    cMat = bone;
        case 'copper';  cMat = copper;
        case 'pink';    cMat = pink;
        case 'lines';   cMat = lines;
        case 'jet';     cMat = jet;
        case 'colorcube';cMat = colorcube;
        case 'viridis'; cMat = viridis;
    end
    warning off
    if nargin == 2
        if nColors>163; error('colorBG function only works for n2 =<63'); end
        cMat = cMat(1:length(cMat)/nColors:length(cMat),:);
    else
        cMat = cMat([1:64],:);
    end

    warning on
end

S = sum(cMat,2);
Duplicates = [];
commonLS = {'-' '--' ':' '-.'};
commonMK = {'none' '+' 's' 'd' '^' 'v'};
N = length(commonMK);
LineStyles ={};Marker ={};
for k = 1:size(cMat,1)
    idx = find(S==S(k));
    if length(idx)==1
        LineStyles{k,1} = commonLS{1};
        Marker{k,1} = commonMK{1};
    else
        if length(idx)/N~= ceil(length(idx)/N)
            idx(end+1:ceil(length(idx)/N)*N)=NaN;
        end
        idx = reshape(idx,N,[]);
        for col = 1:size(idx,2)
            for row = 1:size(idx,1)
                if ~isnan(idx(row,col))
                    LineStyles{idx(row,col),1} = commonLS{col};
                    Marker{idx(row,col),1} = commonMK{row};
                end
            end
        end
    end
    idx = reshape(idx,[],1);idx(isnan(idx))=[];
    S(idx) = NaN;
    Duplicates(idx,1)=k;
end

% --------------------------------------------------------------------------------------------------------------- %
function makeMyFigureNice

set(gcf,'Color',[1 1 1]);
grid off
fig=gcf;
N =  length(fig.Children);
for ii = 1:N
    set(fig.Children(ii),'box', 'off')
    if contains(class(fig.Children(ii)),'Axes')
        fig.Children(ii).FontName = 'Arial';
        fig.Children(ii).Title.FontWeight = 'Normal';
        fig.Children(ii).FontSize = 7;                                                                              % font size

        for ax = 1:length(fig.Children(ii).Children)
            if contains(class(fig.Children(ii).Children(ax)),'matlab.graphics.chart.primitive.Line')
                fig.Children(ii).Children(ax).LineWidth = 2;
            end
        end
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,n,Alpha,Color, MakerSize)
% %% Description - Basilio Goncalves (2020)
% https://www.researchgate.net/profile/Basilio_Goncalves
%
% create a scatter plot and a line of best fit with a shaded area as CI
%
%INPUT
%   x = x axis values
%   y = y axis values
%   n = polynomial degree (Default = 1)
%   Alpha = alpha level betwen 0 and 1(Default = 0.05)
%   Color = color of the plots (Default = black)
%   MakerSize = size of the individual markers (Default = 10)
%-------------------------------------------------------------------------
%OUTPUT
%   rsquared = cell matrix maximum Force value
%   pvalue = cell matrix maximum Force value
%   p1 = handle of the plot

if  nargin < 3 || isempty(n)
    n= length(x);
end
if nargin < 4 || isempty(Alpha)
    Alpha= 0.05;
end

if  nargin < 5 || isempty(Color)
    Color= 'k';
end

if  nargin < 6 || isempty(MarkerSize)
    MakerSize= 10;
end

% delete nan
IDXnan = any(isnan([x y]),2);
x(IDXnan)=[];
y(IDXnan)=[];


% calculate polynomial and line that best fit
[p,S,mu] = polyfit(x,y,n);
[y_fit,delta] = polyval(p,x,S,mu);

% calculate rsquared and p-value
[c, pvalue,rlo,rup] = corrcoef(x,y);
rsquared = c(1,2)^2;
pvalue = pvalue(1,2);
rlo = rlo(1,2);
rup = rup(1,2);

% calculate t-statistic
t = tinv(1-Alpha/2,S.df);

% calculate confidence interval
uB = y_fit+ t*delta;
lB = y_fit- t*delta;


p = sortrows([x y y_fit uB lB],1);                                  %plot data
x= p(:,1)';
y= p(:,2)';
y_fit= p(:,3)';
uB= p(:,4)';
lB= p(:,5)';

hold on
p1 = plot(x,y,'o','MarkerSize',MakerSize,'MarkerFaceColor',Color,'MarkerEdgeColor','none');

p1(end+1) = plot(x,y_fit,'Color',Color,'LineWidth',2);              % plot mean trendline

X=[x,fliplr(x)];                                                    % create continuous x value row vector for plotting
Y=[lB fliplr(uB)];                                                  % create y values for out and then back (Remove NaN)
f1 = fill(X,Y,'r');
alpha 0.2                                                           % transparency
set(f1,'FaceColor', Color,'EdgeColor','none')

% --------------------------------------------------------------------------------------------------------------- %
function plotMuscleForces(Results,Type)
if nargin < 2
    Type ='Normalised';
end

legs = {'left','right'};
for iLeg = 1:2
    leg = legs{iLeg};
    muscle_forces = Results.SO.session1.left;
    muscle_names = fields(muscle_forces);
    
    % only leg side muscles
    muscle_names = muscle_names(endsWith(muscle_names,['_' leg(1)]));
    % remove fileds that contain the term "peaks"
    muscle_names = muscle_names(~contains(muscle_names,'peak_'));                                                   
    line_colors = getColor('viridis',3);

    % remove '_l' or '_r' from the muscle names
    muscle_names_short = {};
    for iMusc = 1:length(muscle_names)
        muscle_names_short{end+1,1} = muscle_names{iMusc}(1:end-2);
        if any(contains(muscle_names_short{end}(end),{'1','2','3'}))                                                % remove numebers 1,2,3
            muscle_names_short{end} = muscle_names_short{end}(1:end-1);
        end
    end
    muscle_names_short = unique(muscle_names_short);

    % create figure with subplots
    n_subplots = length(muscle_names_short);
    [ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);
  
    % loop through all muscles
    for iMusc = 1:length(muscle_names_short)
        indData = {};
        MeanForcesAllSessions = [];
        SDForcesAllSessions = [];

        % assign muscle forces for each session to one 
        for iSess = 1:3
            session = ['session' num2str(iSess)];

            % select either absolute or normalised
            switch Type
                case 'Normalised'
                    muscle_forces = Results.SO_normalised.(session).(leg);
                case 'Absolute'
                    muscle_forces = Results.SO.(session).(leg);
            end

            
            current_muscle = muscle_names_short{iMusc};

            % get all muscle segments for each muscle name
            segments = muscle_names(contains(muscle_names,current_muscle));                                         
            single_muscle_forces = muscle_forces.(segments{1});

            % if a column has all zeros make it all NaN
            single_muscle_forces = ZeroToNaN(single_muscle_forces);
            
            % average the muscle forces for each segment (for each column /trial)
            for iSeg = 2:length(segments)
                single_muscle_forces = (single_muscle_forces + muscle_forces.(segments{iSeg}))./2;                       
            end
            
            if iSess == 3
                single_muscle_forces(:,7)=NaN;
            end


            % plot mean and SD
            indData{iSess} = removeNaNrows(single_muscle_forces,2);
            MeanForcesAllSessions(:,iSess) = nanmean(single_muscle_forces,2);
            SDForcesAllSessions(:,iSess) = nanstd(single_muscle_forces,0,2);

        end

        % run and save SPM plots
        [SPM] = ttest2(indData);
        suptitle([current_muscle ' ' leg])
        
        savedir = [get_main_dir() fp 'SPM_results'];
        if ~isfolder(savedir); mkdir(savedir); end
        saveas(gcf, [savedir fp current_muscle ' ' leg '.jpeg'])
        close(gcf)
        
        % plot force-time curves and add SPM lines on the plot
        axes(ha(iMusc)); hold on
        p = plotShadedSD(MeanForcesAllSessions,SDForcesAllSessions,line_colors);
        
        add_spm_to_plot(SPM)
       
        % change ylim and yticks
        ylim([0 140])
        

        % add ylable to first col
        if any(iMusc == FirstCol)
            yl = ylabel('% Max isom force');
            yl.Rotation = 0;
            yl.HorizontalAlignment ="right";
        end

        % title
        t = title(current_muscle, 'Interpreter','none');
        t.VerticalAlignment = 'top';
    end
    
    % add ticks to the plots
    tight_subplot_ticks (ha,LastRow,FirstCol)                                                                       
    
    % add overall title for the figure 
    suptitle(['muscle forces ' leg])

    % add legend
    lg = legend({'pre' 'sd' 'post' '' 'td' ''});
    lg.FontSize = 12;
    lg.Position = [0.8 0.12 0.05 0.1];

    % make figure nice (backgorund color, font size and type, etc...)
    makeMyFigureNice
    
    % save figure
    saveas(gcf, [savedir fp 'muscle_forces' leg '.jpeg'])

end

% --------------------------------------------------------------------------------------------------------------- %
function p = plotShadedSD(YData,SD,COLOR,Xvalues)
% plot each column on the Mean matrix with each column on the SD matrix

if size(YData,2)~=size(SD,2)
    error('Number of columans in Mean and SD inputs must agree')
end

% [ cMat, cStruct, cNames] = getColorSet(30); % color blind friendly
if nargin<3 || isempty(COLOR)
    cMat = getColor(0,size(YData,2)); % color scheme 2 (Bas)
else
    cMat = COLOR;
end
%line styles
style = {'-','--',':','-.'};
for k = 1:ceil(size(YData,2)/4)
    style = [style style];
end

for ii = 1:size(YData,2)
    hold on
    y1=YData(:,ii)';                             % create main curve
    if nargin <4
        x=1:length(y1);                             % initialize x row vector
    else
        x = Xvalues';
    end

    p(ii) = plot(x,y1,'LineWidth',1);
    set(p(ii),'Color',cMat(ii,:),'LineStyle',style{ii})
    color = p(ii).Color;

    Top= y1+SD(:,ii)';                          % create top of shaed area
    Bottom = y1-SD(:,ii)';                      % create bottom of shaded
    X=[x,fliplr(x)];                            % create continuous x value array for plotting
    Y=[Bottom fliplr(Top)];                     % create y values for out and then back
    f1 = fill(X,Y,color);
    alpha 0.3
    set(f1,'FaceColor', color,'EdgeColor','none')

end

% --------------------------------------------------------------------------------------------------------------- %
function [ha, pos,FirstCol,LastRow,LastCol] = tight_subplot(Nh, Nw, gap, marg_h, marg_w,Size)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos, FirstCol, LastRow] = tight_subplot(Nh, Nw, gap, marg_h, marg_w,Size)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%       Size     4 dimesions of the figure to be plotted.
%                Size == 0 -> make full size figure
%                Size > 0 (eg 0.4) ->  make size figure 0.4 times full size
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplotBG(3,2,[.01 .03],[.1 .01],[.01 .01],[100 200 1200 600])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
%
% Edited by Basilio Goncalves (2021)
% see also: numSubplots

if nargin < 2 || Nw==0
    RemoveSubPlots  = 1;
    LastPlot        = Nh;
    [N,]            = numSubplots(Nh);
    Nh              = N(1);
    Nw              = N(2);
else
    RemoveSubPlots  = 0;
    LastPlot        = Nh*Nw;
end


if nargin<3 || isempty(gap); gap = .04; end
if nargin<4 || isempty(marg_h); marg_h = .04; end
if nargin<5 || isempty(marg_w); marg_w = .04; end
if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end
figure
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh;
% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);

    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);

if nargin<6 || length(Size)==1 && Size == 0 || length(Size)==1 && Size > 1
    set(gcf,'units','normalized','outerposition',[0 0.05 1 0.95])
elseif length(Size)<3 && Size>0
    dimensions = [1*Size 0.95*Size];
    origin = (0 +(1-dimensions))/2;
    set(gcf,'units','normalized','outerposition',[origin dimensions])
elseif length(Size)==4 && ~any(Size>1)
    set(gcf,'units','normalized','outerposition',[Size])
else
    set(gcf, 'Position', Size);
end

FirstCol = 1:Nw:LastPlot; LastRow = Nh*Nw-Nw+1:Nh*Nw; LastCol = Nw:Nw:LastPlot;

if RemoveSubPlots~=0
    delete(ha(LastPlot+1:end));
    ha(LastPlot+1:end) =[];
    LastRow = ceil(LastPlot/Nw);
    LastRow = LastRow*Nw-Nw+1 : LastPlot;
end

% --------------------------------------------------------------------------------------------------------------- %
function tight_subplot_ticks (ha,xt,yt)
% tight_subplot_ticks adds xticklabels and yticklables to the assigned
% subplots in "ha" (if xt OR yt == 1, add labels to all subplots)
%
% tight_subplot_ticks (ha,xt,yt)
% xt = index of subplots to add xticks to
% yt = index of subplots to add yticks to
%
% see also tight_subplotBG

N = length(ha);
for i = 1:N
    axes(ha(i));
    if any(xt==i) || any(xt==0)
        xticklabels(xticks)
    else
        xticklabels('')
    end

    if any(yt==i) || any(yt==0)
        yticklabels(yticks)
    else
        yticklabels('')
    end
end







% --------------------------------------------------------------------------------------------------------------- %
% ------------------------------------------- STATS FUNCTIONS --------------------------------------------------- %
% --------------------------------------------------------------------------------------------------------------- %
function [SPM] = ttest2(indData)

n_groups = length(indData);
combinations = nchoosek([1:n_groups],2); 
Alpha = 0.05;%/size(combinations,1);

ha = tight_subplot(1,n_groups,[],[],[],[0.05 0.3 0.9 0.5]);
SPM = struct;
count = 0;
for iComb = combinations'
    count = count + 1;
    Dataset1 = indData{iComb(1)};
    Dataset2 = indData{iComb(2)};
    
    % run SPM tests
    spmi = spm1d.stats.ttest2(Dataset1',Dataset2');
    spmi = spmi.inference(Alpha);
    
    comparison_name = [num2str(iComb(1)) '_VS_' num2str(iComb(2))];

    % SPM plot with p-values and t-values threshold
    axes(ha(count))
    spmi.plot; spmi.plot_p_values;
    spmi.plot_threshold_label;
    title(comparison_name,'Interpreter','none')

    % from the axis, find the index of the siginifant points and it's p-value
    LinesPlot = ha(count).Children; 
    significant_idx=[]; 
    p_value_idx=[];
    for iLine = 1:length(LinesPlot)

        %  find indexes of patches (signficand shaded areas)   
        if contains(class(LinesPlot(iLine)),'Patch'); significant_idx(end+1) = iLine; end      

        %  find indexes of tesxt with 'P = ' 
        if contains(class(LinesPlot(iLine)),'Text') && contains(LinesPlot(iLine).String,'p '); p_value_idx(end+1) = iLine; end      
    end

    % add spmi results to final struct
    SPM.(['comp_' comparison_name]) = struct;
    flds_spmi = fields(spmi);
    for i = 1:length(flds_spmi)
        curr_fld = flds_spmi{i};
        SPM.(['comp_' comparison_name]).(curr_fld) = spmi.(curr_fld);    
    end

    % add the significance vectors (idx and 
    SPM.(['comp_' comparison_name]).sig_idx = significant_idx;
    SPM.(['comp_' comparison_name]).p_idx = p_value_idx;
end
tight_subplot_ticks(ha,0,0)

% --------------------------------------------------------------------------------------------------------------- %
function add_spm_to_plot(SPM,colors,yPosition)

comparisons = fields(SPM);
nComp = length(comparisons);

% if no colors are selected 
if nargin < 2
    colors = getColor('viridis',nComp);
elseif length(colors) ~= length(comparisons) 
    warnign on 
    warning('number of olors and number of comparions do not match. Auto colors ')
    colors = getColor('viridis',nComp);
end

if nargin < 3
    yPosition = 110;
end

for i = 1:nComp

    % get the x positions of significant differences for each comparisons
    % and create a similar vector y at height 110
    x = SPM.(comparisons{i}).sig_idx;
    x = [1:10, 20:34, 50:68];
    
    % Find the indices where consecutive values change & split the vector into sections
    diff_indices = find(diff(x) ~= 1);
    sections = mat2cell(x, 1, diff([0, diff_indices, numel(x)]));

    % plot each section the sections
    for i = 1:numel(sections)
        x = sections{i};
        y = repmat(110, size(x));

        % Define the position and size of the rectangle
        x_rect = x(1);
        y_rect = yPosition;
        width = length(x);
        height = 0.1;

        % Draw the rectangle
        rectangle('Position', [x_rect, y_rect, width, height], 'FaceColor', colors(i,:), 'EdgeColor','none')
    end
end


% --------------------------------------------------------------------------------------------------------------- %
function x_new = insertDecimalPoints(x,N_decimals)
x_new = [];
if isempty(x)
    return
end

for i = 1:length(x)-1
    x_new = [x_new, x(i), linspace(x(i)+0.00001, x(i+1)-0.00001, N_decimals)];
end
x_new = [x_new, x(end)];







