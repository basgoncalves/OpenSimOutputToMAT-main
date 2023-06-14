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
                
                % create single session results struct
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
    Results = create_muscle_groups(Results);
%     Results = normalise_muscle_forces_max_isom(Results);
        Results = add_participant_mass_to_results_struct(Results);
    Results = normalise_muscle_forces_BW(Results);
    




    % save data
    cd(get_main_dir)
    save('results.mat', 'Results')
end
%% stats
if run_stats

    load([get_main_dir fp 'results.mat'])

    Results = calculatePeaks(Results);

%     plotMuscleForces(Results,'Normalised_BW',1)

    plotContactForces(Results, 1)

%     plotIK(Results,'Absolute', 1)

%     plotID(Results, 1)
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
function Results = add_participant_mass_to_results_struct(Results)
cd(get_main_dir)
P = load('participants_characteristics.mat');

sessions = fields(Results.Mass);

for i = 1:length(sessions)
    current_session = sessions{i};

    % create an all nan row for current session
    n_subjects = length(Results.Subject.(current_session));
    Results.Mass.(current_session) = NaN(1,n_subjects);

    % find the subjects that have data for current sessions
    session_subjects_cols = find(Results.Subject.(current_session));

    for iSubj = 1:length(session_subjects_cols)
        results_col = session_subjects_cols(iSubj);
        Results.Mass.(current_session)(1,results_col) = P.part_char.(current_session).mass(1, iSubj);
    end
end

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
function [out,RemovedIdx] = removeNaNrows (data, Dim)
RemovedIdx = [];
out = data;

if nargin < 2
    Dim = 1;
end

if Dim == 1
    for iRow = 1:size(out,Dim)
        if sum(isnan(out(iRow,:)))==length(out(iRow,:))
            RemovedIdx(end+1) = iCol;
        end

    end
    out(RemovedIdx,:)=[];

elseif Dim == 2
    for iCol = 1:size(out,Dim)
        if sum(isnan(out(:,iCol)))==length(out(:,iCol))
            RemovedIdx(end+1) = iCol;
        end
    end
    out(:,RemovedIdx)=[];

else
    error('Dim input varible must be 1 (rows) or 2 (cols)')
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

        % resultant and peak hip contact forces
        Results.JRL.(session).(leg).hip_fresultant = sum3D(Results.JRL.(session).(leg), ['hip_' leg(1) '_on_pelvis_in_pelvis']);
        Results.JRL.(session).(leg).peak_HCF = calcMaxMovAverage(Results.JRL.(session).(leg).hip_fresultant);

        % knee
        Results.JRL.(session).(leg).knee_fresultant = sum3D(Results.JRL.(session).(leg), ['knee_' leg(1) '_on_tibia_' leg(1) '_in_tibia_' leg(1)]);
        Results.JRL.(session).(leg).peak_KCF = calcMaxMovAverage(Results.JRL.(session).(leg).knee_fresultant);

        % ankle
        Results.JRL.(session).(leg).ankle_fresultant = sum3D(Results.JRL.(session).(leg), ['ankle_' leg(1) '_on_talus_' leg(1) '_in_talus_' leg(1)]);
        Results.JRL.(session).(leg).peak_ACF = calcMaxMovAverage(Results.JRL.(session).(leg).ankle_fresultant);
        
        % peak muscle forces
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
function A = ZeroToNaN(A)

% Loop through each column
for col = 1:size(A, 2)
    if all(A(:, col) == 0)

        % Replace column with NaN values
        A(:, col) = NaN;
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function muscle_names_unique = find_unique_names_no_numbers(muscle_names)

% Remove the last two elements from each cell using cellfun and indexing
muscle_names_unique = cellfun(@(x) x(1:end-2), muscle_names, 'UniformOutput', false);

% remove numbers from the end
for i = 1:numel(muscle_names_unique)
    iMusc = muscle_names_unique{i};
    if ~isempty(str2num(iMusc(end)))
        muscle_names_unique{i} = muscle_names_unique{i}(1:end-1);
    end
end
muscle_names_unique = unique(muscle_names_unique);

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
function plotMuscleForces(Results,Type, Bonf_corr) % plot all muscle forces and run t-test SPM to compare (pre,post, & TD)
if nargin < 2
    Type ='Normalised';
end

legs = {'left','right'};
line_colors = getColor('parula',3);

% muscle names
muscle_names = fields(Results.SO.session1.left);
muscle_names = muscle_names(endsWith(muscle_names,'_l')); % only leg side muscles
muscle_names = muscle_names(~contains(muscle_names,'peak_')); % remove fileds that contain the term "peaks"

% muscle_names_unique = find_unique_names_no_numbers(muscle_names);
muscle_groups = {'glut_max'; 'glut_med'; 'glut_min';'adductors'; 'hamstrings'; 'iliopsoas'; 'tfl'; 'ext_rot'; 'rect_fem'}; 

% create figure with subplots
% n_subplots = length(muscle_names_unique);
n_subplots = length(muscle_groups);
[ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);

mean_differences = struct;
% loop through all muscles
for iMusc = 1:length(muscle_groups)
%     for iMusc = 1:length(muscle_names_unique)
    indData = {[] [] []};
    MeanForcesAllSessions = [];
    SDForcesAllSessions = [];
    % assign muscle forces for each session to one
    for iSess = 1:3
        session = ['session' num2str(iSess)];

        for iLeg = 1:2
            leg = legs{iLeg};

            muscle_names = cellfun(@(str) [str(1:end-1) leg(1)], muscle_names, 'UniformOutput', false);

            % select either absolute or normalised
            switch Type
                case 'Normalised_max_iso_force'
                    muscle_forces = Results.SO_normalised.(session).(leg);
                case 'Normalised_BW'
                    muscle_forces = Results.SO_normalised_BW.(session).(leg);
                case 'Absolute'
                    muscle_forces = Results.SO.(session).(leg);
            end

            current_muscle = [muscle_groups{iMusc} '_' leg(1)];
            %             current_muscle = [current_muscle '_' leg(1)];

            % get all muscle segments for each muscle name
            %             segments = muscle_names(contains(muscle_names,current_muscle));
            %             single_muscle_forces = muscle_forces.(segments{1});
            single_muscle_forces = muscle_forces.(current_muscle);

            % if a column has all zeros make it all NaN
            single_muscle_forces = ZeroToNaN(single_muscle_forces);

            % average the muscle forces for each segment (for each column /trial)
%             for iSeg = 2:length(segments)
%                 single_muscle_forces = (single_muscle_forces + muscle_forces.(segments{iSeg}))./2;
%             end

%             if iSess == 3
%                 single_muscle_forces(:,7)=NaN;
%             end

            % remove columns with just NaNs
            single_muscle_forces = removeNaNrows(single_muscle_forces,2);

            % find dimensions of the data
            [nRows,nCols] = size(single_muscle_forces);
            idxCols = [size(indData{iSess},2)+1 : size(indData{iSess},2)+nCols];

            % add individual participant data to a cell
            indData{iSess}(:,idxCols) = single_muscle_forces;

%             figure()
%             plot(single_muscle_forces)
%             title([session ': ' current_muscle ' ' segments{1} leg],'Interpreter', 'None')

        end

        % calculate mean and SD for each group
        MeanForcesAllSessions(:,iSess) = nanmean(indData{iSess},2);
        SDForcesAllSessions(:,iSess) = nanstd(indData{iSess},0,2);
    end
    
    % run and save SPM plots
    group_types = [1,1,2];
    [SPM] = ttest_spm(indData,group_types);
    suptitle([current_muscle])

    savedir = [get_main_dir() fp 'SPM_results'];
    if ~isfolder(savedir); mkdir(savedir); end
    saveas(gcf, [savedir fp current_muscle ' ' leg '.jpeg'])
    close(gcf)

    % plot force-time curves and add SPM lines on the plot
    axes(ha(iMusc)); hold on
    p = plotShadedSD(MeanForcesAllSessions,SDForcesAllSessions,line_colors);

    % change ylim and yticks
    %     ylim([0 180])
    xlim([0 101])

    r = add_spm_to_plot(SPM,140);

    short_muscle = current_muscle(1:end-2);
    mean_differences.(short_muscle) = {};
    differences_loc.(short_muscle) = {};
    comparisons = fields(SPM);
    nComparisons = length(comparisons);
    for iComp = 1:nComparisons
        current_comp = comparisons{iComp};
        mean_differences.(short_muscle){iComp,1} = current_comp;
        mean_differences.(short_muscle){iComp,2} = SPM.(current_comp).MeanDifference;
                differences_loc.(short_muscle){iComp,1} = current_comp;
        differences_loc.(short_muscle){iComp,2} = SPM.(current_comp).sig_idx;
    end

    % add ylable to first col
    if any(iMusc == FirstCol)
        yl = ylabel('muscle forces [N/BW]');
        yl.Rotation = 0;
        yl.HorizontalAlignment ="right";
    end

    % title
    t = title(current_muscle(1:end-2), 'Interpreter','none');
    t.VerticalAlignment = 'top';
    t.FontName = 'Arial';
    t.FontSize = 10;
end

% add ticks to the plots
if contains(Type,'Normalised')
    tight_subplot_ticks (ha,LastRow,FirstCol)
else
    tight_subplot_ticks (ha,LastRow,0)
end

% add overall title for the figure
% suptitle(['muscle forces ' Type])

% add legend
lg = legend({'pre' 'sd' 'post' '' 'td' '' 'pre vs post' 'pre vs TD' 'post vs TD'});
lg.FontSize = 10;
lg.Location = 'best';

% make figure nice (backgorund color, font size and type, etc...)
makeMyFigureNice

% save figure
saveas(gcf, [savedir fp 'muscle_forces_BW' leg '.jpeg'])
save('mean_differences_m.mat', 'mean_differences')
save('differences_loc_m.mat', 'differences_loc')


% --------------------------------------------------------------------------------------------------------------- %
function plotContactForces(Results, Bonf_corr) % plot all muscle forces and run t-test SPM to compare (pre,post, & TD)

legs = {'left','right'};
line_colors = getColor('parula',3);

% JRF names containing fx,fy,fz,fresultant
JRL_names = fields(Results.JRL.session1.left);
JRL_names = [JRL_names(endsWith(JRL_names,{'_fx','_fy','fz','_fresultant'}))]; 

% replace "_l" and "_r_" with "_leg_" to plot just the variable for a
% single leg
JRL_names = strrep(JRL_names,'_r_','_leg_');
JRL_names = strrep(JRL_names,'_l_','_leg_');
JRL_names = unique(JRL_names);

% create figure with subplots
n_subplots = length(JRL_names);
[ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);

% loop through all muscles
for iJRL = 1:length(JRL_names)
    indData = {[] [] []};
    MeanForcesAllSessions = [];
    SDForcesAllSessions = [];
    % assign muscle forces for each session to one
    for iSess = 1:3
        session = ['session' num2str(iSess)];
     
        for iLeg = 1:2
            leg = legs{iLeg};

            % select either absolute or normalised
            contact_forces = Results.JRL.(session).(leg);
            current_JRL_variable =  strrep(JRL_names{iJRL},'_leg_',['_' leg(1) '_']);

            % get all muscle segments for each muscle name
            single_varibale_JRL = contact_forces.(current_JRL_variable);

            % if a column has all zeros make it all NaN
            single_varibale_JRL = ZeroToNaN(single_varibale_JRL);

            % get mass per session (do not move outside the loop to avoid
            % deleting first participant more than once)
            mass_current_session = Results.Mass.(session)(~isnan(Results.Mass.(session)));

            % participant not good data
%             if iSess == 3
%                 single_varibale_JRL(:,7)=NaN;
%                 mass_current_session(1) = [];
%                 warning on
%                 warning ('deleting data participant in column 7 for session 3')
%             end

            % remove columns with just NaNs
            single_varibale_JRL = removeNaNrows(single_varibale_JRL,2);

            % flip the contact forces for the left leg
            if contains(leg(1),'l') && contains(current_JRL_variable(end-1:end),'fz')
                single_varibale_JRL = -single_varibale_JRL;
            end
            
            % normalise to BW
            single_varibale_JRL_BW_normalised = single_varibale_JRL./(mass_current_session*9.81);

            % find dimensions of the data
            [nRows,nCols] = size(single_varibale_JRL_BW_normalised);
            idxCols = [size(indData{iSess},2)+1 : size(indData{iSess},2)+nCols];

            % add individual participant data to a cell
            indData{iSess}(:,idxCols) = single_varibale_JRL_BW_normalised;

        end

        % calculate mean and SD for each group
        MeanForcesAllSessions(:,iSess) = nanmean(indData{iSess},2);
        SDForcesAllSessions(:,iSess) = nanstd(indData{iSess},0,2);
    end

    parts_JRL_variable = split(current_JRL_variable,'_');
    current_JRL_variable_title = [parts_JRL_variable{1} '_' parts_JRL_variable{end}];

    % run and save SPM plots
    group_types = [1,1,2];
    [SPM] = ttest_spm(indData,group_types);
    suptitle([current_JRL_variable_title])

    savedir = [get_main_dir() fp 'SPM_results' fp 'JRL'];
    if ~isfolder(savedir); mkdir(savedir); end
    saveas(gcf, [savedir fp current_JRL_variable_title '.jpeg'])
    close(gcf)

    % plot force-time curves and add SPM lines on the plot
    axes(ha(iJRL)); hold on
    p = plotShadedSD(MeanForcesAllSessions,SDForcesAllSessions,line_colors);

    % change ylim and yticks
%     ylim([0 6])
    
    r = add_spm_to_plot(SPM,[],8.5);

% short_muscle = current_muscle(1:end-2);
    mean_differences.(JRL_names{iJRL}) = {};
    differences_loc.(JRL_names{iJRL}) = {};
    comparisons = fields(SPM);
    nComparisons = length(comparisons);
    for iComp = 1:nComparisons
        current_comp = comparisons{iComp};
        mean_differences.(JRL_names{iJRL}){iComp,1} = current_comp;
        mean_differences.(JRL_names{iJRL}){iComp,2} = SPM.(current_comp).MeanDifference;
        differences_loc.(JRL_names{iJRL}){iComp} = SPM.(current_comp).sig_idx;
    end


    % add ylable to first col
    if any(iJRL == FirstCol)
        yl = ylabel('Joint contact forces (BW)');
        yl.Rotation = 0;
        yl.HorizontalAlignment ="right";
    end

    % title
    t = title(current_JRL_variable_title, 'Interpreter','none');
    t.VerticalAlignment = 'top';
end

% add ticks to the plots
tight_subplot_ticks (ha,LastRow,0)

% add legend
lg = legend({'pre' 'sd' 'post' '' 'td' '' 'pre vs post' 'pre vs TD' 'post vs TD'});
lg.FontSize = 10;
lg.FontName = 'Arial';
% lg.Position = [0.8 0.12 0.05 0.1];
lg.Location = 'best';

% make figure nice (backgorund color, font size and type, etc...)
makeMyFigureNice

% save figure
saveas(gcf, [savedir fp 'JointReactionLoads.jpeg'])
save('mean_differences_jrf.mat', 'mean_differences')
save('differences_loc_jrf.mat', 'differences_loc')



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
function [SPM] = ttest_spm(indData,group_types)
% type: 1 = independent / 2 = paired


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
    if group_types(iComb(1)) ==  group_types(iComb(2))
        spmi = spm1d.stats.ttest_paired(Dataset1',Dataset2'); % if both datasets are from the same group
        type = 'paired';
    else
        spmi = spm1d.stats.ttest2(Dataset1',Dataset2'); % if datasets are from different groups
        type = '2-sample';
    end

    spmi = spmi.inference(Alpha);

    comparison_name = [num2str(iComb(1)) 'VS' num2str(iComb(2))];

    % SPM plot with p-values and t-values threshold
    axes(ha(count))
    spmi.plot; spmi.plot_p_values;
    spmi.plot_threshold_label;
    title([comparison_name ' ' type],'Interpreter','none')

    % from the axis, find the index of the siginifant points and it's p-value
    LinesPlot = ha(count).Children;
    significant_idx={};
    p_value_idx=[];
    MeanDifference ={};
    for iLine = 1:length(LinesPlot)

        %  find indexes of patches (signficand shaded areas)
        if contains(class(LinesPlot(iLine)),'Patch')
            xData = unique(round(LinesPlot(iLine).XData));
            cols = [length(significant_idx)+1:length(significant_idx) + length(xData)];
            significant_idx{1,end+1} = xData'; 

            % mean and CI for the area of each patch
            [MD,lCI,uCI] = meanDif (Dataset1(xData,:),Dataset2(xData,:),Alpha,type,0); 
            
            % create a string type "MD (lCI-uCI)" with the correct decimal points
            MeanDifference{1,end+1} = [determine_decimal_points(MD) ' (' determine_decimal_points(lCI) '-' determine_decimal_points(uCI) ')'];
        end

        %  find indexes of tesxt with 'P = '
        if contains(class(LinesPlot(iLine)),'Text') && contains(LinesPlot(iLine).String,'p ')
            p_value_idx(end+1) = str2num(LinesPlot(iLine).String(4:end));
        end
    end

    % add spmi results to final struct
    SPM.(['comp_' comparison_name]) = struct;
    flds_spmi = fields(spmi);
    for i = 1:length(flds_spmi)
        curr_fld = flds_spmi{i};
        SPM.(['comp_' comparison_name]).(curr_fld) = spmi.(curr_fld);
    end

    % flip variables because the first elements represent last patches on
    % the figure 
    p_value_idx = flip(p_value_idx);
    significant_idx = flip(significant_idx);
    MeanDifference = flip(MeanDifference);

    % add the significance vectors (idx and
    SPM.(['comp_' comparison_name]).sig_idx = significant_idx;
    SPM.(['comp_' comparison_name]).p_idx = p_value_idx;
    SPM.(['comp_' comparison_name]).MeanDifference = MeanDifference;
end
tight_subplot_ticks(ha,0,0)

% --------------------------------------------------------------------------------------------------------------- %
function r = add_spm_to_plot(SPM,colors,yPosition)
r = [];
comparisons = fields(SPM);
nComp = length(comparisons);

% if no colors ddre selected
if nargin < 2 || isempty(colors)
    colors = getColor('parula',nComp);
elseif length(colors) ~= length(comparisons)
    warning on
    warning('number of olors and number of comparions do not match. Auto colors ')
    colors = getColor('parula',nComp);
end
colors_rect = getColor('gray',3);

% y coordinates of the sig rectangle
if nargin < 3
    yPosition = max(ylim) * 1.2;
    ylim([min(ylim) ceil(yPosition)])
end

distance_between_rectanges = yPosition*0.05;

for iComp = 1:nComp

    % get the x positions of significant differences for each comparisons
    % indx_significant_differences = [1:50,60:85]; % test values
    sections_significant_differences = SPM.(comparisons{iComp}).sig_idx;

    % Find the indices where consecutive values change & split the vector into sections
    if isempty(sections_significant_differences)
        disp('no sig dif found')
        continue
    end

    % plot each section the sections
    for i = 1:numel(sections_significant_differences)
        indx_significant_differences = sections_significant_differences{i};

        % Define the position and size of the rectangle
        x_rect = indx_significant_differences(1);
        y_rect = yPosition-distance_between_rectanges*iComp;
        width = length(indx_significant_differences);
        height = distance_between_rectanges/2;

        % Draw the rectangle
        if ~isnan(x_rect)
            r = rectangle('Position', [x_rect, y_rect, width, height], 'FaceColor', colors_rect(iComp,:), 'EdgeColor','none');
        end

    end

    % Draw fake rectangle for the legend entry
    hold on
    rectangle_plot = plot(NaN, NaN, 's', 'MarkerFaceColor', colors_rect(iComp,:), 'MarkerEdgeColor', 'none');
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


% --------------------------------------------------------------------------------------------------------------- %
function plotIK(Results,Type, Bonf_corr) % plot kinematics and run t-test SPM to compare (pre,post, & TD)
if nargin < 2
    Type ='Absolute';
end

legs = {'left','right'};
line_colors = getColor('parula',3);

% angle names
angle_names = fields(Results.IK.session1.left);
angle_names = angle_names(endsWith(angle_names,'_l')); % only leg side angles
angle_names = angle_names(~contains(angle_names,{'mtp','subtalar'})); % remove fields that contain the term "peaks"

angle_names_unique = find_unique_names_no_numbers(angle_names);

% create figure with subplots
n_subplots = length(angle_names_unique);
[ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);

% loop through all angles
for iAng= 1:length(angle_names_unique)
    indData = {[] [] []};
    MeanForcesAllSessions = [];
    SDForcesAllSessions = [];
    % assign angle for each session to one
    for iSess = 1:3
        session = ['session' num2str(iSess)];

        for iLeg = 1:2
            leg = legs{iLeg};

            % select either absolute or normalised
            switch Type
                case 'Normalised'
                    angle_forces = Results.IK_normalised.(session).(leg);
                case 'Absolute'
                    angle_forces = Results.IK.(session).(leg);
            end

            current_angle = angle_names_unique{iAng};
            current_angle = [current_angle '_' leg(1)];


            % get all angle segments for each angle name
%             segments = angle_names(contains(angle_names,current_angle));
            single_angle_forces = angle_forces.(current_angle);

            % if a column has all zeros make it all NaN
            single_angle_forces = ZeroToNaN(single_angle_forces);

            % average the angle forces for each segment (for each column /trial)
%             for iSeg = 2:length(segments)
%                 single_angle_forces = (single_angle_forces + angle_forces.(segments{iSeg}))./2;
%             end

%             if iSess == 3
%                 single_angle_forces(:,7)=NaN;
%             end

            % remove columns with just NaNs
            single_angle_forces = removeNaNrows(single_angle_forces,2);

            % find dimensions of the data
            [nRows,nCols] = size(single_angle_forces);
            idxCols = [size(indData{iSess},2)+1 : size(indData{iSess},2)+nCols];

%             figure()
%             plot(single_angle_forces)
%             title([session ': ' current_angle],'Interpreter', 'None')

            % add individual participant data to a cell
            indData{iSess}(:,idxCols) = single_angle_forces;

        end

        % calculate mean and SD for each group
        MeanAnglesAllSessions(:,iSess) = nanmean(indData{iSess},2);
        SDAnglesAllSessions(:,iSess) = nanstd(indData{iSess},0,2);
    end
    % run and save SPM plots
    group_types = [1,1,2];
    [SPM] = ttest_spm(indData,group_types, Bonf_corr);
    suptitle([current_angle])

    savedir = [get_main_dir() fp 'SPM_results'];
    if ~isfolder(savedir); mkdir(savedir); end
    saveas(gcf, [savedir fp current_angle ' ' leg '.jpeg'])
    close(gcf)

    % plot force-time curves and add SPM lines on the plot
    axes(ha(iAng)); hold on
    p = plotShadedSD(MeanAnglesAllSessions,SDAnglesAllSessions,line_colors);

    % change ylim and yticks
%     ylim([0 140])
    
    r = add_spm_to_plot(SPM,140);

    % add ylable to first col
    if any(iAng == FirstCol)
        yl = ylabel('Angle [°]');
        yl.Rotation = 0;
        yl.HorizontalAlignment ="right";
    end

    % title
    t = title(current_angle, 'Interpreter','none');
    t.VerticalAlignment = 'top';
end

% add ticks to the plots
if contains(Type,'Normalised')
    tight_subplot_ticks (ha,LastRow,FirstCol)
else
    tight_subplot_ticks (ha,LastRow,0)
end

% add overall title for the figure
suptitle(['kinematics ' Type])

% add legend
lg = legend({'pre' 'sd' 'post' '' 'td' '' 'pre vs post' 'pre vs TD' 'post vs TD'});
lg.FontSize = 12;
lg.Position = [0.8 0.12 0.05 0.1];

% make figure nice (backgorund color, font size and type, etc...)
makeMyFigureNice

% save figure
saveas(gcf, [savedir fp 'kinematics' leg '.jpeg'])


% --------------------------------------------------------------------------------------------------------------- %
function plotID(Results, Bonf_corr) % plot all moments and run t-test SPM to compare (pre,post, & TD)

legs = {'left','right'};
line_colors = getColor('parula',3);

% JRF names containing fx,fy,fz,fresultant
ID_names = fields(Results.ID.session1.left);
ID_names = ID_names(contains(ID_names,{'moment'})); % remove fields that contain the term "peaks"
ID_names = ID_names(~contains(ID_names,{'pelvis','lumbar','mtp','subtalar'})); % remove fields that contain the term "peaks"

% JRL_names = [JRL_names(endsWith(JRL_names,{'_fx','_fy','fz','_fresultant'}))]; 

% replace "_l" and "_r_" with "_leg_" to plot just the variable for a
% single leg
ID_names = strrep(ID_names,'_r_','_leg_');
ID_names = strrep(ID_names,'_l_','_leg_');
ID_names = unique(ID_names);

% create figure with subplots
n_subplots = length(ID_names);
[ha, ~,FirstCol,LastRow,~] = tight_subplot(n_subplots,0,[0.008 0.01],[0.05 0.02],[0.08 0.01],0.99);

% loop through all moments
for iID = 1:length(ID_names)
    indData = {[] [] []};
    MeanMomentsAllSessions = [];
    SDMomentsAllSessions = [];
    % assign muscle forces for each session to one
    for iSess = 1:3
        session = ['session' num2str(iSess)];
     
        for iLeg = 1:2
            leg = legs{iLeg};

            % select either absolute or normalised
            contact_forces = Results.ID.(session).(leg);
            current_ID_variable =  strrep(ID_names{iID},'_leg_',['_' leg(1) '_']);

            % get all muscle segments for each muscle name
            single_variable_ID = contact_forces.(current_ID_variable);

            % if a column has all zeros make it all NaN
            single_variable_ID = ZeroToNaN(single_variable_ID);

            % get mass per session (do not move outside the loop to avoid
            % deleting first participant more than once)
            mass_current_session = Results.Mass.(session)(~isnan(Results.Mass.(session)));

            % participant not good data
%             if iSess == 3
%                 single_varibale_ID(:,7)=NaN;
%                 mass_current_session(1) = [];
%                 warning on
%                 warning ('deleting data participant in column 7 for session 3')
%             end

            % remove columns with just NaNs
            single_variable_ID = removeNaNrows(single_variable_ID,2);

            % flip the contact forces for the left leg
%             if contains(leg(1),'l') && contains(current_ID_variable(end-1:end),'fz')
%                 single_varibale_ID = -single_varibale_ID;
%             end
            
            % normalise to BW
            single_varibale_ID_BW_normalised = single_variable_ID./(mass_current_session*9.81);

            % find dimensions of the data
            [nRows,nCols] = size(single_varibale_ID_BW_normalised);
            idxCols = [size(indData{iSess},2)+1 : size(indData{iSess},2)+nCols];

            % add individual participant data to a cell
            indData{iSess}(:,idxCols) = single_varibale_ID_BW_normalised;

        end

        % calculate mean and SD for each group
        MeanMomentsAllSessions(:,iSess) = nanmean(indData{iSess},2);
        SDMomentsAllSessions(:,iSess) = nanstd(indData{iSess},0,2);
    end

    parts_ID_variable = split(current_ID_variable,'_r_');
    current_ID_variable_title = [parts_ID_variable{1} '_' parts_ID_variable{end}];
    parts_ID_variable = split(current_ID_variable_title,'_');
    current_ID_variable_title = [parts_ID_variable{1} ' ' parts_ID_variable{2} ' ' parts_ID_variable{end}];

    % run and save SPM plots
    group_types = [1,1,2];
    [SPM] = ttest_spm(indData,group_types, Bonf_corr);
    suptitle([current_ID_variable_title])

    savedir = [get_main_dir() fp 'SPM_results' fp 'ID'];
    if ~isfolder(savedir); mkdir(savedir); end
    saveas(gcf, [savedir fp current_ID_variable_title '.jpeg'])
    close(gcf)

    % plot force-time curves and add SPM lines on the plot
    axes(ha(iID)); hold on
    p = plotShadedSD(MeanMomentsAllSessions,SDMomentsAllSessions,line_colors);

    % change ylim and yticks
    ylim([-0.15 0.15])
    xlim([0 100])

    spm_colors = getColor('hot',3);
    
    r = add_spm_to_plot(SPM,spm_colors,0.1);

    % add ylable to first col
    if any(iID == FirstCol)
        yl = ylabel('Joint Moments [N/BW]');
        yl.Rotation = 0;
        yl.HorizontalAlignment ="right";
    end

    % title
    t = title(current_ID_variable_title, 'Interpreter','none');
    t.VerticalAlignment = 'top';
end

% add ticks to the plots
tight_subplot_ticks (ha,LastRow,0)

% add legend
lg = legend({'pre' 'sd' '' 'post' '' 'td' '' 'pre vs post' 'pre vs TD' 'post vs TD'});
lg.FontSize = 12;
lg.Position = [0.8 0.12 0.05 0.1];

% make figure nice (backgorund color, font size and type, etc...)
makeMyFigureNice

% save figure
saveas(gcf, [savedir fp 'JointMoments.jpeg'])

% --------------------------------------------------------------------------------------------------------------- %
function Results = normalise_muscle_forces_BW(Results)

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
        muscle_names = fields(Results.SO.session1.right);
        muscle_names = muscle_names(~contains(muscle_names,{'reserve', 'calcn','time', 'duration'}));
        muscle_names_r = muscle_names(endsWith(muscle_names,'_r'));
        muscle_names_l = cellfun(@(str) [str(1:end-1) 'l'], muscle_names_r, 'UniformOutput', false);

        for iMusc = 1:length(muscle_names_r)
%             current_muscle = muscle_names{iMusc};

            % try to add
            try
                right_norm_temp.(muscle_names_r{iMusc})(:,iCol) = Results.SO.(session).right.(muscle_names_r{iMusc})(:,iCol) / Results.Mass.(session)(iCol)*9.81  ;
                left_norm_temp.(muscle_names_l{iMusc})(:,iCol) = Results.SO.(session).left.(muscle_names_l{iMusc})(:,iCol) / Results.Mass.(session)(iCol)*9.81;
            catch
                right_norm_temp.(muscle_names_r{iMusc}) = [];
                left_norm_temp.(muscle_names_l{iMusc}) = [];
                right_norm_temp.(muscle_names_r{iMusc})(:,iCol) = Results.SO.(session).right.(muscle_names_r{iMusc})(:,iCol) / Results.Mass.(session)(iCol)*9.81;
                left_norm_temp.(muscle_names_l{iMusc})(:,iCol) = Results.SO.(session).left.(muscle_names_l{iMusc})(:,iCol) / Results.Mass.(session)(iCol)*9.81;
            end

        end
    end

    % append data to the main Results struct
    Results.SO_normalised_BW.(session).right = right_norm_temp;
    Results.SO_normalised_BW.(session).left = left_norm_temp;
end


% ---------------------------------------------------------------------- %
% -------------------------------------------------------------------------%
function Results_BW = create_muscle_groups(Results_BW)

legs = {'left'; 'right'};
    sessions = {'session1'; 'session2'; 'session3'};
%     joints = {'HCF'; 'KCF'; 'ACF'};

 musc_seg = {'glut_max'; 'glut_med'; 'glut_min'; 'add_mag'}; 

 for i = 1:size(sessions,1)
     for j = 1:size(legs,1)
         for k = 1:size(musc_seg,1)
             Results_BW.SO.(sessions{i}).(legs{j}).([musc_seg{k} '_' legs{j}(1)]) = Results_BW.SO.(sessions{i}).(legs{j}).([musc_seg{k} '1_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).([musc_seg{k} '2_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).([musc_seg{k} '3_' legs{j}(1)]);
%              Results_BW.SO.(sessions{i}).(legs{j}).(['peak_' musc_seg{k} '_' legs{j}(1) '_val']) = Results_BW.SO.(sessions{i}).(legs{j}).(['peak_' musc_seg{k} '1_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_' musc_seg{k} '2_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_' musc_seg{k} '3_' legs{j}(1) '_val']);

         end
         Results_BW.SO.(sessions{i}).(legs{j}).(['hamstrings_' legs{j}(1)]) = Results_BW.SO.(sessions{i}).(legs{j}).(['semimem_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['semiten_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['bifemlh_' legs{j}(1)]);
         Results_BW.SO.(sessions{i}).(legs{j}).(['adductors_' legs{j}(1)]) = Results_BW.SO.(sessions{i}).(legs{j}).(['add_mag_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['add_long_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['add_brev_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['pect_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['grac_' legs{j}(1)]);
         Results_BW.SO.(sessions{i}).(legs{j}).(['iliopsoas_' legs{j}(1)]) = Results_BW.SO.(sessions{i}).(legs{j}).(['iliacus_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['psoas_' legs{j}(1)]);
         Results_BW.SO.(sessions{i}).(legs{j}).(['ext_rot_' legs{j}(1)]) = Results_BW.SO.(sessions{i}).(legs{j}).(['gem_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['quad_fem_' legs{j}(1)]) + Results_BW.SO.(sessions{i}).(legs{j}).(['peri_' legs{j}(1)]);

%          Results_BW.SO.(sessions{i}).(legs{j}).(['peak_hamstrings_' legs{j}(1) '_val']) = Results_BW.SO.(sessions{i}).(legs{j}).(['peak_semimem_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_semiten_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_bifemlh_' legs{j}(1) '_val']);
%          Results_BW.SO.(sessions{i}).(legs{j}).(['peak_adductors_' legs{j}(1) '_val']) = Results_BW.SO.(sessions{i}).(legs{j}).(['peak_add_mag_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_add_long_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_add_brev_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_pect_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_grac_' legs{j}(1) '_val']);
%          Results_BW.SO.(sessions{i}).(legs{j}).(['peak_iliopsoas_' legs{j}(1) '_val']) = Results_BW.SO.(sessions{i}).(legs{j}).(['peak_iliacus_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_psoas_' legs{j}(1) '_val']);
%          Results_BW.SO.(sessions{i}).(legs{j}).(['peak_ext_rot_' legs{j}(1) '_val']) = Results_BW.SO.(sessions{i}).(legs{j}).(['peak_gem_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_quad_fem_' legs{j}(1) '_val']) + Results_BW.SO.(sessions{i}).(legs{j}).(['peak_peri_' legs{j}(1) '_val']);
     end
 end

% gluteus maximus
% gluteus medius
% gluteus minimus
% hamstrings (semimembranosus, semitendinosus, and biceps femoris long head)
% adductors (adductor magnus, longus, and brevis, pectineus and gracilis)
% iliopsoas
% rectus femoris
% tensor fasciae latae
% deep external rotators (gemellus superior, gemellus inferior, quadratus femoris and piriformis)

% -------------------------------------------------------------------------------%

function [MD,LB,UB] = meanDif (D1,D2, Alpha,Type,ConvertToPercentage)
% calc mean difference in percentage
%Type: 1(default) = paired; 2 = independent;
% ConvertToPercentage: 1(default) = true; 0 = false;
if ~exist('Alpha')||isempty(Alpha);  Alpha = 0.05; end

if ~exist('Type')||isempty(Type)
    Type = 1; 
elseif contains(Type,'2-sample') || contains(Type,'independent')
    Type = 2;

elseif contains(Type,'1-sample') || contains(Type,'paired')
    Type = 1;
end

if ~exist('ConvertToPercentage')||isempty(ConvertToPercentage)
    ConvertToPercentage = 1; 
end

mean1 = mean(D1,1);
mean2 = mean(D2,1);

n1 = size (mean1,2);
n2 = size (mean2,2);

if Type == 1
    if ConvertToPercentage == 1
        D_diff = (mean2 - mean1)./mean1*100;
    elseif ConvertToPercentage == 0
        D_diff = mean2 - mean1;
    end
    
    MD = mean(D_diff);
    df = n1-1;
    CIdiff = std(D_diff)/sqrt(n1)*tinv(1-Alpha/2,df);
    LB = MD-CIdiff;
    UB = MD+CIdiff;
    
elseif Type == 2      % https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals5.html
    
    MD = mean(mean2)-mean(mean1);
    
    std1 = std(mean1);
    std2 = std(mean2);
    
    df = n1+n2-2;
    t = tinv(1-Alpha/2,df);
    var1= (n1-1)*std1^2;
    var2= (n2-1)*std2^2;
    pooledSE = sqrt((var1 + var2) / df) * sqrt(1/n1+1/n2);
    
    CIdiff = pooledSE*tinv(1-Alpha/2,df);
    LB = MD-CIdiff;
    UB = MD+CIdiff;
    
end

% -------------------------------------------------------------------%

function [formattedValue,decimalPoints] = determine_decimal_points(value)

if nargin < 1
    value = 34;
end

magnitude = abs(value);
if magnitude >= 1000
    decimalPoints = 0;  % No decimal points for large values
elseif magnitude >= 10
    decimalPoints = 1;  % One decimal point for medium values
else
    decimalPoints = 2;  % Two decimal points for small values
end

formattedValue = sprintf(['%.' num2str(decimalPoints) 'f'], value);




