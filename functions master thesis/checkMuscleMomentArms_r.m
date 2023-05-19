function momentArmsAreWrong = checkMuscleMomentArms_r(modelFilename, motionFilename)
%CHECKMUSCLEMOMENTARM Summary of this function goes here
%   Detailed explanation goes here
tic

import org.opensim.modeling.*

motion = Storage(motionFilename);
model = Model(modelFilename);

model.initSystem();
state = model.initSystem();

flexIndexR = model.getCoordinateSet().getIndex('hip_flexion_r');
flexCoordR = model.updCoordinateSet().get(flexIndexR);
rotIndexR = model.getCoordinateSet().getIndex('hip_rotation_r');
rotCoordR = model.updCoordinateSet().get(rotIndexR);
addIndexR = model.getCoordinateSet().getIndex('hip_adduction_r');
addCoordR = model.updCoordinateSet().get(addIndexR);
flexIndexRknee = model.getCoordinateSet().getIndex('knee_angle_r');
flexCoordRknee = model.updCoordinateSet().get(flexIndexRknee);
flexIndexRank = model.getCoordinateSet().getIndex('ankle_angle_r');
flexCoordRank = model.updCoordinateSet().get(flexIndexRank);

numMuscles = model.getMuscles().getSize();
muscleIndices = []; muscleNames = {};
for i = 0 : numMuscles - 1
    tmp_muscleName = char(model.getMuscles().get(i).getName());
    if contains(tmp_muscleName, 'add') || contains(tmp_muscleName, 'gl') ...
            || contains(tmp_muscleName, 'semi') || contains(tmp_muscleName, 'bf') ...
            || contains(tmp_muscleName, 'grac') || contains(tmp_muscleName, 'piri') ...
            || contains(tmp_muscleName, 'sart') || contains(tmp_muscleName, 'tfl') ...
            || contains(tmp_muscleName, 'iliacus') || contains(tmp_muscleName, 'psoas') ...
            || contains(tmp_muscleName, 'rect')
        if contains(tmp_muscleName,'_r')

            muscleIndices = [muscleIndices, i];
            muscleNames{end+1} = tmp_muscleName;
        end
    end
end

flexMomentArms = zeros(motion.getSize(), length(muscleIndices));
addMomentArms = zeros(motion.getSize(), length(muscleIndices));
rotMomentArms = zeros(motion.getSize(), length(muscleIndices));

for i = 1:motion.getSize()
    flexAngleR = motion.getStateVector(i-1).getData().get(flexIndexR) / 180 * pi;
    rotAngleR = motion.getStateVector(i-1).getData().get(rotIndexR) / 180 * pi;
    addAngleR = motion.getStateVector(i-1).getData().get(addIndexR) / 180 * pi;
    flexAngleRknee = motion.getStateVector(i-1).getData().get(flexIndexRknee) / 180 * pi;
    flexAngleRank = motion.getStateVector(i-1).getData().get(flexIndexRank) / 180 * pi;

    % Update the state with the joint angle
    coordSet = model.updCoordinateSet();
    coordSet.get(flexIndexR).setValue(state, flexAngleR);
    coordSet.get(rotIndexR).setValue(state, rotAngleR);
    coordSet.get(addIndexR).setValue(state, addAngleR);
    coordSet.get(flexIndexRknee).setValue(state, flexAngleRknee);
    coordSet.get(flexIndexRank).setValue(state, flexAngleRank);

    % Realize the state to compute dependent quantities
    model.computeStateVariableDerivatives(state);
    model.realizeVelocity(state);

    % Compute the moment arm for each muscle
    for j = 1:length(muscleIndices)
        muscleIndex = muscleIndices(j);
        if strcmp(muscleNames{j}(end), 'r')
            flexMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, flexCoordR);
            flexMomentArms(i,j) = flexMomentArm;
            rotMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, rotCoordR);
            rotMomentArms(i,j) = rotMomentArm;
            addMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, addCoordR);
            addMomentArms(i,j) = addMomentArm;
            flexMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, flexCoordRknee);
            flexMomentArms(i,j) = flexMomentArm;
            rotMomentArm = model.getMuscles().get(muscleIndex).computeMomentArm(state, flexCoordRank);
            rotMomentArms(i,j) = rotMomentArm;
        end
    end
end

toc

threshold = 0.005;
discontinuity = [];
fDistC = figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
hold on;
legArr = {};
for i = 1 : size(flexMomentArms, 2)
    dy = diff(flexMomentArms(:, i));
    discontinuity_indices = find(abs(dy) > threshold);
    if size(discontinuity_indices, 1) > 0
        fprintf(2, ['Discontinuity detected at ' muscleNames{i} ' at flexion moment arm \n'])
        plot(flexMomentArms(:, i))
        plot(discontinuity_indices, flexMomentArms(discontinuity_indices, i), 'rx');
        discontinuity = [discontinuity, i];
        legArr = [legArr; [muscleNames{i} ' flexion']; ' '];
    end

    dy = diff(addMomentArms(:, i));
    discontinuity_indices = find(abs(dy) > threshold);
    if size(discontinuity_indices, 1) > 0
        fprintf(2, ['Discontinuity detected at ' muscleNames{i} ' at adduction moment arm \n'])
        plot(addMomentArms(:, i))
        plot(discontinuity_indices, addMomentArms(discontinuity_indices, i), 'rx');
        discontinuity = [discontinuity, i];
        legArr = [legArr; [muscleNames{i} ' adduction']; ' '];
    end

    dy = diff(rotMomentArms(:, i));
    discontinuity_indices = find(abs(dy) > threshold);
    if size(discontinuity_indices, 1) > 0
        fprintf(2, ['Discontinuity detected at ' muscleNames{i} ' at rotation moment arm \n'])
        plot(rotMomentArms(:, i))
        plot(discontinuity_indices, rotMomentArms(discontinuity_indices, i), 'rx');
        discontinuity = [discontinuity, i];
        legArr = [legArr; [muscleNames{i} ' rotation']; ' '];
    end
end

title(motionFilename);

if size(discontinuity, 1) > 0
    legend(legArr);
    
    ylabel('Muscle Moment Arms with discontinuities (m)');
    xlabel('Frame (after start time)');
    drawnow;
    fprintf(2, '\n\nYou should alter the model - most probably you have to reduce the radius of corresponding wrap objects for the identified muscles\n\n\n')
    momentArmsAreWrong = 1;
else
    close(fDistC);
    disp('No discontinuities detected');
    momentArmsAreWrong = 0;
end

figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
plot(flexMomentArms);
title(['All muscle moment arms in motion ' motionFilename]);
legend(muscleNames, 'Interpreter', 'none');
ylabel('Hip Flexion Moment Arm (m)');
xlabel('Frame (after start time)');
drawnow;
end
