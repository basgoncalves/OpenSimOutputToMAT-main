function [frameZero,cycle,tempData, frequency] = get_frame_zero(stoFile)

currentFolder = fileparts(fileparts(fileparts(stoFile)));
load(fullfile(currentFolder, 'settings.mat'));

tempData = load_sto_file(stoFile);
% s.cycle.right.start =  s.cycle.right.start + 3;
% calculate the first event
if tempData.time(1) >0
    beginning_recording_new = round(tempData.time(1)/(1/frequency));
    preframes = preframes - (beginning_recording_new-1);
    if isfield(cycle, 'left') && isfield(cycle, 'right')
        cycle.left.start = cycle.left.start - (beginning_recording_new-1);
        cycle.left.end = cycle.left.end - (beginning_recording_new-1);
        cycle.left.footOff = cycle.left.footOff - (beginning_recording_new-1);
        cycle.right.start = cycle.right.start - (beginning_recording_new-1);
        cycle.right.end = cycle.right.end - (beginning_recording_new-1);
        cycle.right.footOff = cycle.right.footOff - (beginning_recording_new-1);
    elseif isfield(cycle, 'left')
        cycle.left.start = cycle.left.start - (beginning_recording_new-1);
        cycle.left.end = cycle.left.end - (beginning_recording_new-1);
        cycle.left.footOff = cycle.left.footOff - (beginning_recording_new-1);
    else
        cycle.right.start = cycle.right.start (beginning_recording_new-1);
        cycle.right.end = cycle.right.end - (beginning_recording_new-1);
        cycle.right.footOff = cycle.right.footOff - (beginning_recording_new-1);
    end
end


if isfield(cycle, 'left') && isfield(cycle, 'right')
    frameZero = min(min(cycle.left.start), min(cycle.right.start)) - 1;
elseif isfield(cycle, 'left')
    frameZero = min(cycle.left.start) - 1;
else
    frameZero = min(cycle.right.start) - 1;
end


if exist('preframes', 'var')
    %                     frameZero = frameZero - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
    frameZero = floor(preframes) - frameZero; % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit
end