function currentFolder = find_c3d_starting_at_zero
paths = {'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\pre\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P01\post\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P02\pre\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P02\post\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P03\pre\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P03\post\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P04\pre\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P04\post\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P05\pre\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\P05\post\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD01\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD04\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD06\output automization';
    'C:\Users\Balu\Nextcloud\Documents\MA\Daten\TD07\output automization';};

for i = 1:length(paths)
    outputPath = paths{i};
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
            if firstFrame > preframes
                return
            end

        end
    end
end

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

