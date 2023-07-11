function saveDataToStruct_bodyKin(outputPath)
% clear;
% outputPath = 'C:\Users\Willi\Documents\UniDataLokal\TorsionToolComparison\Simulationsresults';

addpath 'C:\Users\Balu\Nextcloud\Documents\MA\Code'
addpath 'C:\Users\Balu\Nextcloud\Documents\MA\Code\MSKmodelling'
activate_msk_modelling
modelList = GetSubDirsFirstLevelOnly(outputPath);

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

ErrorScores = zeros(4,1);
%% 
for p = 1 : numel(modelList)
    disp(' ');
    probandFolder = [outputPath filesep modelList{p}];
    trialList = GetSubDirsFirstLevelOnly(probandFolder);
    for trialNr = 1 : numel(trialList)
        disp(['currently processing:' ' ' modelList{p} ' ' '|' ' ' trialList{trialNr}]);
        
        currentFolder = [probandFolder filesep trialList{trialNr}];

        clear preframes EMG EMGChannels ratioEmgToFrames emgFrequency;
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
        if cycleIsSorted
            
            ikFile = fullfile(currentFolder, 'Output', 'IK', 'IK.mot');
            if isfile(ikFile)


%                 fid = fopen(fullfile(currentFolder, 'Output', 'IK', 'ikSettings.xml'), 'r');
%                 f=fread(fid,'*char')';
%                 fclose(fid);
%                 i0 = strfind(f, '<marker_file>');
%                 i1 = strfind(f, '</marker_file>');
                
%                 markerFileName = f(i0 + 13 : i1 - 1);
                [~,parentFolder] = fileparts(currentFolder);
               if contains(parentFolder, 'Dynamic22') || contains(parentFolder, 'Dynamic09')
                   continue
               end
                markerFileName = [fileparts(fileparts(fileparts(currentFolder))) '\C3D\' parentFolder '\marker_experimental.trc'];

%                 markerFileName = strrep(markerFileName, '/Users/niklasnagel/Documents/Study/GaitAnalyseData', 'C:\Users\Willi\Documents\UniDataLokal\TorsionToolComparison\GaitAnalyseData');
                
                markerData = load_marker_trc(markerFileName);
                markerNames = fieldnames(markerData);
                heelMarker = markerNames(contains(markerNames, 'HEE'));
                heelMarker = heelMarker{1}(1 : end-2);
                heelMarker = strrep(heelMarker, 'LH', 'H');
                heelMarker = strrep(heelMarker, 'RH', 'H');
                
                leftHeel = [cell2mat(markerData.(['L' heelMarker '_X'])), cell2mat(markerData.(['L' heelMarker '_Y'])), cell2mat(markerData.(['L' heelMarker '_Z']))];
                rightHeel = [cell2mat(markerData.(['R' heelMarker '_X'])), cell2mat(markerData.(['R' heelMarker '_Y'])), cell2mat(markerData.(['R' heelMarker '_Z']))];


                %crop data from left and right leg to stance phase
                clear tempData tempFieldNames tempStructLeft tempStructRight frameZero;
                tempData = load_sto_file(ikFile);
                tempFieldNames = fieldnames(tempData);
                tempStructLeft = struct;
                tempStructRight = struct;
                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end
                
                if exist('preframes', 'var')
                    preframes = preframes + 1;
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        clear tempStructLeft;
                        for i = 1 : numel(tempFieldNames)
                            tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                        end
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).footOffFrame = double(cycle.left.footOff(j) - cycle.left.start(j));

                        heelPositionStrike1 = rightHeel(cycle.left.start(j)+1, :);
                        heelPositionStrike2 = rightHeel(cycle.left.end(j), :);
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).walkingDirection = heelPositionStrike2 - heelPositionStrike1;
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        clear tempStructRight;
                        for i = 1 : numel(tempFieldNames)
                            tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j)+1 - frameZero : cycle.right.end(j) - frameZero - 1);
                        end
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)+1) / frequency;
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).footOffFrame = double(cycle.right.footOff(j) - cycle.right.start(j)+1);

                        heelPositionStrike1 = rightHeel(cycle.right.start(j)+1, :);
                        heelPositionStrike2 = rightHeel(cycle.right.end(j), :);
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).walkingDirection = heelPositionStrike2 - heelPositionStrike1;
                    end
                end
                
                % check for IK_fehler.txt bzw IK_withTorsoOutOfRecommendation.txt
                % ErrorScore 1: Alles cool
                % ErrorScore 2: Nur Marker errors am Oberkörper zu hoch
                % ErrorScore 3: Marker error bei Inverse kinematic zu hoch
                % ErrorScore 4: Fehler bei Simulation
                
                % if isfile([currentFolder '/IK_fehler.txt']) && isfile([currentFolder '/IK_withTorsoOutOfRecommendation.txt']) && ~isfile([currentFolder '/IK_noProb.txt'])
                %     % beide files vorhanden --> ErrorScore 3 --> bad
                %     ErrorScores(3,1) = ErrorScores(3,1)+1;
                %     data.ErrorScoreList.ErrorScore3{ErrorScores(3,1)} = currentFolder;
                if  isfile([currentFolder '/IK_withTorsoOutOfRecommendation.txt']) && ~isfile([currentFolder '/IK_fehler.txt'])  %elseif
                    % nur IK_withTorso... vorhanden --> ErrorScore 2 --> ok
                    ErrorScores(2,1) = ErrorScores(2,1)+1;
                    data.ErrorScoreList.ErrorScore2{ErrorScores(2,1)} = currentFolder;
                else
                    % nichts vorhanden --> ErrorScore 1 --> perfect
                    ErrorScores(1,1) = ErrorScores(1,1)+1;
                    data.ErrorScoreList.ErrorScore1{ErrorScores(1,1)} = currentFolder;
                end
            end
            
            
            idFile = fullfile(currentFolder, 'Output', 'ID', 'inverse_dynamics.sto');
            if isfile(idFile)
                data.ID.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                clear tempData tempFieldNames tempStructLeft tempStructRight frameZero;
                tempData = load_sto_file(idFile);
                tempFieldNames = fieldnames(tempData);
                tempStructLeft = struct;
                tempStructRight = struct;
                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end
                
                if exist('preframes', 'var')
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        clear tempStructLeft;
                        for i = 1 : numel(tempFieldNames)
                            tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                        end
                        data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                        data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        clear tempStructRight;
                        for i = 1 : numel(tempFieldNames)
                            tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                        end
                        data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                        data.ID.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                    end
                end
            end
            
            
            soFile = fullfile(currentFolder, 'Output', 'SO', '_StaticOptimization_force.sto');
            if isfile(soFile)
                data.SO.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                clear tempData tempFieldNames tempStructLeft tempStructRight frameZero;
                tempData = load_sto_file(soFile);
                tempFieldNames = fieldnames(tempData);
                tempStructLeft = struct;
                tempStructRight = struct;
                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end
                
                
                if exist('preframes', 'var')
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        clear tempStructLeft;
                        for i = 1 : numel(tempFieldNames)
                            tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                        end
                        data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                        data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        clear tempStructRight;
                        for i = 1 : numel(tempFieldNames)
                            tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                        end
                        data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                        data.SO.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                    end
                end
            end
            
            soFile = fullfile(currentFolder, 'Output', 'SO', '_StaticOptimization_activation.sto');
            if isfile(soFile)
                data.SO_Activation.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                clear tempData tempFieldNames tempStructLeft tempStructRight frameZero;
                tempData = load_sto_file(soFile);
                tempFieldNames = fieldnames(tempData);
                tempStructLeft = struct;
                tempStructRight = struct;
                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end
                
                
                if exist('preframes', 'var')
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        clear tempStructLeft;
                        for i = 1 : numel(tempFieldNames)
                            tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                        end
                        data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                        data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        clear tempStructRight;
                        for i = 1 : numel(tempFieldNames)
                            tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                        end
                        data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                        data.SO_Activation.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                    end
                end
            end
            
            jrlFile = fullfile(currentFolder, 'Output', 'JRL', '_JointReaction_ReactionLoads.sto');
            if isfile(jrlFile)
                data.JRL.(modelList{p}).model_mass = model_mass;
                %crop data from left and right leg to stance phase
                clear tempData tempFieldNames tempStructLeft tempStructRight frameZero;
                tempData = load_sto_file(jrlFile);
                tempFieldNames = fieldnames(tempData);
                tempStructLeft = struct;
                tempStructRight = struct;
                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end
                
                
                if exist('preframes', 'var')
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end
                
                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        clear tempStructLeft;
                        for i = 1 : numel(tempFieldNames)
                            tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                        end
                        data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = tempStructLeft;
                        data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).durationInSeconds = double(cycle.left.end(j) - cycle.left.start(j)) / frequency;
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        clear tempStructRight;
                        for i = 1 : numel(tempFieldNames)
                            tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                        end
                        data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = tempStructRight;
                        data.JRL.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).durationInSeconds = double(cycle.right.end(j) - cycle.right.start(j)) / frequency;
                    end
                end
            elseif isfile(ikFile) % jrl file doesn't exist --> error in simulation --> add to ErrorScore 4
                ErrorScores(4,1) = ErrorScores(4,1)+1;
                data.ErrorScoreList.ErrorScore4{ErrorScores(4,1)} = currentFolder;
            end



            bodyKinFile = fullfile(currentFolder, 'Output', 'BodyKin', '_BodyKinematics_pos_global.sto');
            if isfile(bodyKinFile)
                %crop data from left and right leg to stance phase
                clear tempData tempFieldNames tempStructLeft tempStructRight frameZero;
                tempData = load_sto_file(bodyKinFile);
                tempFieldNames = fieldnames(tempData);
                tempStructLeft = struct;
                tempStructRight = struct;
                if isfield(cycle, 'left') && isfield(cycle, 'right')
                    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
                elseif isfield(cycle, 'left')
                    frameZero = min(cycle.left.start) - 1;
                else
                    frameZero = min(cycle.right.start) - 1;
                end


                if exist('preframes', 'var')
                    frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
                end

                if isfield(cycle, 'left')
                    for j = 1 : size(cycle.left.start, 2)
                        clear tempStructLeft;
                        for i = 1 : numel(tempFieldNames)
                            tempStructLeft.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.left.start(j) - frameZero : cycle.left.end(j) - frameZero - 1);
                        end

                        orientYfieldnames = tempFieldNames(contains(tempFieldNames, 'Oy'));
                        walkingDir = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']).walkingDirection;
                        walkingDirAngle = atan2d(walkingDir(3), walkingDir(1));
                        for i = 1 : numel(orientYfieldnames)
                            tempStructLeft.(orientYfieldnames{i}) = tempStructLeft.(orientYfieldnames{i}) - walkingDirAngle;
                        end

                        tempStructLeft = rmfield(tempStructLeft, 'time');
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']) = mergestructs(data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_left']), tempStructLeft);
                    end
                end
                if isfield(cycle, 'right')
                    for j = 1 : size(cycle.right.start, 2)
                        clear tempStructRight;
                        for i = 1 : numel(tempFieldNames)
                            tempStructRight.(tempFieldNames{i}) = tempData.(tempFieldNames{i})(cycle.right.start(j) - frameZero : cycle.right.end(j) - frameZero - 1);
                        end

                        orientYfieldnames = tempFieldNames(contains(tempFieldNames, 'Oy'));
                        walkingDir = data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']).walkingDirection;
                        walkingDirAngle = atan2d(walkingDir(3), walkingDir(1));
                        for i = 1 : numel(orientYfieldnames)
                            tempStructRight.(orientYfieldnames{i}) = tempStructRight.(orientYfieldnames{i}) - walkingDirAngle;
                        end

                        tempStructRight = rmfield(tempStructRight, 'time');
                        data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']) = mergestructs(data.IK.(modelList{p}).(['T_' trialList{trialNr} '_' num2str(j) '_right']), tempStructRight);
                    end
                end
            end
       
        else
            disp('reprocess this trial! cycle variable is unsorted!');
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
delete([outputPath '/bodyKinematics.mat']);
dataFile = [outputPath '/bodyKinematics.mat'];
save(dataFile, 'data', 'Scores_table');