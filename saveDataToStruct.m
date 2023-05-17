%%
function saveDataToStruct(outputPath)
if nargin < 1
    outputPath = 'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P04\post\output automization';
end
modelList = GetSubDirsFirstLevelOnly(outputPath);

ErrorScores = zeros(4,1);
%%
for p = 1 : numel(modelList)
    disp(' ');
    probandFolder = [outputPath filesep modelList{p}];
    trialList = GetSubDirsFirstLevelOnly(probandFolder);
    for trialNr = 1 : numel(trialList)
        disp(['currently processing:' ' ' modelList{p} ' ' '|' ' ' trialList{trialNr}]);
        currentFolder = [probandFolder filesep trialList{trialNr}];

        if ~isCycleSorted(currentFolder)
            disp('reprocess this trial! cycle variable is unsorted!');
            continue
        end
        load(fullfile(currentFolder, 'settings.mat'));

        %% IK
        ikFile = fullfile(currentFolder, 'Output', 'IK', 'IK.mot');
        if isfile(ikFile)
            %crop data from left and right leg to stance phase
            [tempData,tempFieldNames,tempStructLeft,tempStructRight,frameZero] = resetTempData(ikFile);

            if isfield(cycle, 'left')
                for j = 1 : size(cycle.left.start, 2)
                    clear tempStructLeft;
                    for i = 1 : numel(tempFieldNames)
                        tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) : cycle.left.end(j) - 1);
                    end
                    data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                    data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                    data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).footOffFrame = double(cycle.left.footOff(j) - cycle.left.start(j));
                end
            end
        end
    end
end

%% create table with #trials per Error Score & percentage of total #trials
RowDescription = {'number of trials', 'percentage'};
totalNumberOfTrials = [sum(ErrorScores); 100];
ErrorScore1(1,1) = ErrorScores(1,1);
ErrorScore2(1,1) = ErrorScores(2,1);
ErrorScore3(1,1) = ErrorScores(3,1);
ErrorScore4(1,1) = ErrorScores(4,1);
ErrorScore1(2,1) = ErrorScores(1,1)/totalNumberOfTrials(1,1) * 100;
ErrorScore2(2,1) = ErrorScores(2,1)/totalNumberOfTrials(1,1) * 100;
ErrorScore3(2,1) = ErrorScores(3,1)/totalNumberOfTrials(1,1) * 100;
ErrorScore4(2,1) = ErrorScores(4,1)/totalNumberOfTrials(1,1) * 100;
Scores_table = table(totalNumberOfTrials, ErrorScore1, ErrorScore2, ErrorScore3, ErrorScore4);

%% overwrite existing .mat file
delete([outputPath '\dataStruct_IK_results.mat']);
dataFile = [outputPath '\dataStruct_IK_results.mat'];
save(dataFile, 'data', 'Scores_table');

%%
function cycleIsSorted = isCycleSorted(currentFolder)

load(fullfile(currentFolder, 'settings.mat'));
cycleIsSorted = 0;
if isfield(cycle, 'left') && isfield(cycle, 'right')
    if issorted(cycle.left.start) && issorted(cycle.right.start)
        cycleIsSorted = 1;
    end
elseif isfield(cycle, 'left')
    if issorted(cycle.left.start)
        cycleIsSorted = 1;
    end
else
    if issorted(cycle.right.start)
        cycleIsSorted = 1;
    end
end

%%
function  [tempData,tempFieldNames,tempStructLeft,tempStructRight,frameZero] = resetTempData(filePath)
tempData = load_sto_file(filePath);
tempFieldNames = fieldnames(tempData);
tempStructLeft = struct;
tempStructRight = struct;
frameZero = get_frame_zero(filePath);

