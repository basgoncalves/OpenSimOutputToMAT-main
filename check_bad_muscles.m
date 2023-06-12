M = load('C:\Git\Kira_MSc_data\P03\pre\output automization\dataStruct_ErrorScores_no_trials_removed.mat');

model_path = 'C:\Git\Kira_MSc_data\P03\pre\Model\FINAL_PERSONALISEDTORSIONS_scaled_final.osim';
import org.opensim.modeling.*

% Load the OpenSim model
model = Model(model_path);

% Get the number of muscles in the model
numMuscles = model.getMuscles().getSize();

% Create an empty cell array to store the muscle names
ankleMuscleNames = cell(0);

% Loop through all the muscles in the model
for i = 0:numMuscles-1
    muscle = model.getMuscles().get(i);
    if contains(char(muscle.getName()),'soleus')
        break
    end

    % Check if the muscle crosses the ankle joint
    if muscle.getGeometryPath().getPathPointSet().contains('ankle')
        ankleMuscleNames{end+1} = char(muscle.getName());
    end
end

% Display the ankle-spanning muscle names
disp(ankleMuscleNames)

% 
% figure;plot(M.data.SO.FINAL_PERSONALISEDTORSIONS_scaled_final.T_Dynamic06_1_left.quad_fem_r)

import org.opensim.modeling.*

% Load the OpenSim model
model = Model('path/to/your/model.osim');

% Get the number of muscles in the model
numMuscles = model.getMuscles().getSize();

% Create a cell array to store the joint names
jointNames = cell(0);

% Loop through all the muscles in the model
for i = 0:numMuscles-1
    muscle = model.getMuscles().get(i);
    
    % Get the muscle's path points
    pathPoints = muscle.getGeometryPath().getPathPointSet();
    
    % Loop through all the path points
    numPathPoints = pathPoints.getSize();
    for j = 0:numPathPoints-1
        pathPoint = pathPoints.get(j);
        
        % Get the joint name associated with the path point
        jointName = char(pathPoint.getParentFrame().getName());
        
        % Check if the joint name is not already in the list
        if ~ismember(jointName, jointNames)
            jointNames{end+1} = jointName;
        end
    end
end

% Display the joint names
disp(jointNames);

