
function stats_Kira()


add_repos_to_path
[Results] = get_data_struct; 
S = get_subjects;

gather_data_in_struct = false;
run_stats = true;

%% organise data in struct 
if gather_data_in_struct
    Analyses = {'IK','ID','SO','SO_Activation','JRL'};
    for iSubj = 1:length(S.subjects)
        for iSess = 1:length(S.sessions)

            Paths = get_data_paths(iSubj,iSess);
            if ~isfolder(Paths.session)
                disp([Paths.session ' does not exist'])
                continue
            end

            load(Paths.matResults)

            for iVar = 1:length(Analyses)
                curr_analysis = Analyses{iVar};
                data_struct = data.(curr_analysis).FINAL_PERSONALISEDTORSIONS_scaled_final;
                results_struct = Results.(curr_analysis).(['session' num2str(iSess)]);

                [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct);

                Results.(curr_analysis).(['session' num2str(iSess)]) = results_struct;

            end
        end
    end
    cd(get_main_dir)
    save('results.mat', 'Results')
end

%% stats
if run_stats
    load([get_main_dir fp 'results.mat'])
    
    Results = calculate_peaks(Results);
    sessions = {'session1'; 'session2'; 'session3'};
    joint = 'KCF';
    sessions = [sessions; sessions];
%     titles = ['Session 1: right HCF';'Session 1: right KCF'; 'Session 2: right HCF'; 'Session 2: right KCF';'Session 3: right HCF'; 'Session 3: right KCF' ];
    figure
    [ax, pos,FirstCol,LastRow,LastCol] = tight_subplotBG(6,0);
    for iDOF = 1:6
    axes(ax(iDOF)); hold on; grid on;
    plot(Results.JRL.(sessions{iDOF}).left.(joint)); hold on;
    plot(Results.JRL.(sessions{iDOF}).left.(['peak_' joint '_loc']), Results.JRL.(sessions{iDOF}).left.(['peak_' joint '_val']),'*')
    title([sessions{iDOF} ' left ' joint])
    if contains(sessions{iDOF},'session3')
        joint = 'HCF';
    end
    if any(iDOF == FirstCol)
        ylabel('Joint contact force [N]')
    elseif any(iDOF == LastRow)
        xlabel('Gait Cycle [%]')
    end
    end
    load_participant_characteristics()

%     x = Results.JRL.session1.left.peak_HCF';
%     y = Results.JRL.session1.right.peak_HCF';
%     [rsquared,pvalue, p1,rlo,rup] = plotCorr (x,y,1,0.05);  

   
  

end

%% -------------------------------------------------------------------------------------------------------------- %
% ---------------------------------------------------- FUCNTIONS ------------------------------------------------ %
% --------------------------------------------------------------------------------------------------------------- %
function out  = fp
out  = filesep;

% --------------------------------------------------------------------------------------------------------------- %
function add_repos_to_path
activeFile = [mfilename('fullpath') '.m'];

OpenSimOutputToMAT_Dir = fileparts(activeFile);
addpath(genpath(OpenSimOutputToMAT_Dir))
disp([OpenSimOutputToMAT_Dir ' added to path'])


msk_modellingDir  = [fileparts(OpenSimOutputToMAT_Dir) '\MSKmodelling'];
rmpath(genpath(msk_modellingDir))
disp([msk_modellingDir ' added to path'])

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
function load_participant_characteristics()
cd(get_main_dir)
load('participants_characteristics.mat')

% --------------------------------------------------------------------------------------------------------------- %
function [Results] = get_data_struct()

Results = struct;
leg = {'right','left'};
for i = 1:3
    for l = 1:numel(leg)
        session = ['session' num2str(i)];
        Results.IK.(session).(leg{l}) = struct;
        Results.ID.(session).(leg{l}) = struct;
        Results.SO.(session).(leg{l}) = struct;
        Results.SO_Activation.(session).(leg{l}) = struct;
        Results.JRL.(session).(leg{l}) = struct;
        Results.Mass.(session) = [];
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function Paths = get_data_paths(iSubj,iSess)

S = get_subjects;
session_path = [get_main_dir fp S.subjects{iSubj} fp S.sessions{iSess} fp 'output automization\FINAL_PERSONALISEDTORSIONS_scaled_final'];
trialNames      = {getfolders(session_path).name};

Paths.main      = get_main_dir;
Paths.session   = session_path;
Paths.trials    = cellfun(@(c)[session_path fp c],trialNames,'uni',false);
Paths.trialNames = trialNames;
Paths.matResults = [fileparts(session_path) fp 'dataStruct_ErrorScores_no_trials_removed.mat'];

% --------------------------------------------------------------------------------------------------------------- %
function [raw,tnorm,results_struct] = time_norm_per_session(data_struct,results_struct)
% data_struct = a struct similar to the output of load_sto_file


% add mass to results_struct AND remove "

trialNames = fields(data_struct);
if contains(trialNames{1},'mass')
    trialNames(1) = [];
end

% remove variables that contain calc
variables  = fields(data_struct.(trialNames{1}));
variables = variables(~contains(variables,'calc')); 

% create structs
raw   = struct; raw.left = struct; raw.right = struct;
tnorm = struct; tnorm.left = struct; tnorm.right = struct;

% crete the full struct;
for iVar = variables'
    raw.right.(iVar{1}) = [];
    raw.left.(iVar{1}) = [];
    tnorm.right.(iVar{1}) = [];
    tnorm.left.(iVar{1}) = [];
end


% separate right and left trials into different structs
for iTrial = trialNames'
    for iVar = variables'
        try
            currentTrialData = data_struct.(iTrial{1}).(iVar{1});  % if variable doesn't exist in one trial, make it NaN
        catch
            currentTrialData = nan(101,1);
        end
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
    if ~isfield(results_struct.right,iVar{1})
        results_struct.right.(iVar{1}) = [];
        results_struct.left.(iVar{1}) = [];
    end
    
    % right leg (try in case the "iVar" is empty)
    try
        results_struct.right.(iVar{1})(:,end+1) = tnorm.right.(iVar{1});
    end
    % leg leg
    try
        results_struct.left.(iVar{1})(:,end+1) = tnorm.left.(iVar{1});
    end
end

% --------------------------------------------------------------------------------------------------------------- %
function Results = calculate_peaks(Results)


Legs = {'right','left'};

for i = 1:3
    for l = 1:length(Legs)
        session = ['session' num2str(i)];
        leg = Legs{l};

        % find peak HCF
        Results.JRL.(session).(leg).HCF = sum3D(Results.JRL.(session).(leg), ['hip_' leg(1) '_on_pelvis_in_pelvis']);

        lim1 = round(length(Results.JRL.(session).(leg).HCF)*0.08);
        lim2 = round(length(Results.JRL.(session).(leg).HCF)*0.35);
        lim3 = round(length(Results.JRL.(session).(leg).HCF)*0.8);

        [Results.JRL.(session).(leg).peak_HCF_val(1,:), Results.JRL.(session).(leg).peak_HCF_loc(1,:)] = calc_max(Results.JRL.(session).(leg).HCF(lim1:lim2,:));
        [Results.JRL.(session).(leg).peak_HCF_val(2,:), Results.JRL.(session).(leg).peak_HCF_loc(2,:)] = calc_max(Results.JRL.(session).(leg).HCF(lim2:lim3,:));
        Results.JRL.(session).(leg).peak_HCF_loc(1,:) = Results.JRL.(session).(leg).peak_HCF_loc(1,:) + lim1-1;
        Results.JRL.(session).(leg).peak_HCF_loc(2,:) = Results.JRL.(session).(leg).peak_HCF_loc(2,:) + lim2-1;

        % find peak KCF
        Results.JRL.(session).(leg).KCF = sum3D(Results.JRL.(session).(leg), ['knee_' leg(1) '_on_tibia_' leg(1) '_in_tibia_' leg(1)]);

        lim1 = round(length(Results.JRL.(session).(leg).KCF)*0.08);
        lim2 = round(length(Results.JRL.(session).(leg).KCF)*0.3);
        lim3 = round(length(Results.JRL.(session).(leg).KCF)*0.8);

        [Results.JRL.(session).(leg).peak_KCF_val(1,:), Results.JRL.(session).(leg).peak_KCF_loc(1,:)] = calc_max(Results.JRL.(session).(leg).KCF(lim1:lim2,:));
        [Results.JRL.(session).(leg).peak_KCF_val(2,:), Results.JRL.(session).(leg).peak_KCF_loc(2,:)] = calc_max(Results.JRL.(session).(leg).KCF(lim2:lim3,:));
        Results.JRL.(session).(leg).peak_KCF_loc(1,:) = Results.JRL.(session).(leg).peak_KCF_loc(1,:) + lim1-1;
        Results.JRL.(session).(leg).peak_KCF_loc(2,:) = Results.JRL.(session).(leg).peak_KCF_loc(2,:) + lim2-1;

        Muscles = fields(Results.SO.(session).(leg));
        for m = 1:length(Muscles)
            muscle = Muscles{m};
            Results.SO.(session).(leg).(['peak_' muscle]) = calc_max(Results.SO.(session).(leg).(muscle));
        end

    end
end

% --------------------------------------------------------------------------------------------------------------- %
function out = sum3D(Results_substruct, variable)
x = Results_substruct.([variable '_fx']);
y = Results_substruct.([variable '_fy']);
z = Results_substruct.([variable '_fz']);

out = sqrt(x.^2 + y.^2 + z.^2);

%------------------ calc max moving average 
function [val,loc] = calc_max(time_curve)

[val, loc] = max(movmean(time_curve,3));

% --------------------------------------------------------------------------------------------------------------- %
function folders = getfolders(directory, contianing_string, IgnoreCase)

folders = dir (directory);
try 
    folders(1:2) =[];                                           % delete "../" and "./"
catch 
    return
end
folders = folders([folders.isdir]);                         % select only rows that contain folders 

if nargin > 1
    if IgnoreCase == 0
        folders = folders(contains({folders.name},contianing_string));
    else
        folders = folders(contains({folders.name},contianing_string,"IgnoreCase",true));
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
function TimeNormalizedData = TimeNorm (Data,fs)

TimeNormalizedData=[];

for col = 1: size (Data,2)
    
    currentData = Data(:,col);
    currentData(isnan(currentData))=[];
    if length(currentData)<3
        TimeNormalizedData(1:101,col)= NaN;
        continue
    end
     
    timeTrial = 0:1/fs:size(currentData,1)/fs;
    timeTrial(end)=[];
    Tnorm = timeTrial(end)/101:timeTrial(end)/101:timeTrial(end);
    
    TimeNormalizedData(1:101,col)= interp1(timeTrial,currentData,Tnorm)';
end

% --------------------------------------------------------------------------------------------- %
function [ha, pos,FirstCol,LastRow,LastCol] = tight_subplotBG(Nh, Nw, gap, marg_h, marg_w,Size)
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

% ------------------------------------------------------------------------- %
function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots. 
%
% Inputs
% n - the desired number of subplots.     
%  
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.     
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p = 
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
%
%
% Rob Campbell - January 2010
   
    
while isprime(n) & n>4, 
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end    
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design 
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
