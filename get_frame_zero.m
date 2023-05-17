function frame_first_c3dEvent_in_openSim_frame = get_frame_zero(stoFile)

currentFolder = fileparts(fileparts(fileparts(stoFile)));
s = load(fullfile(currentFolder, 'settings.mat'));

tempData = load_sto_file(stoFile);
% s.cycle.right.start =  s.cycle.right.start + 3;
% calculate the first event
if isfield(s.cycle, 'left') && isfield(s.cycle, 'right')
    frame_first_c3dEvent_in_openSim_frame = min(min(s.cycle.left.start), min(s.cycle.right.start)) - 1;
elseif isfield(s.cycle, 'left')
    frame_first_c3dEvent_in_openSim_frame = min(s.cycle.left.start) - 1;
else
    frame_first_c3dEvent_in_openSim_frame = min(s.cycle.right.start) - 1;
end

frame_firstEventMokka = frame_first_c3dEvent_in_openSim_frame + s.firstFrame + 2; % plus 2 = 1 for the frameZero calculation and 1 when calculating the czcle.left/right.start   
difference_btw_first_frame_and_first_event = frame_first_c3dEvent_in_openSim_frame - s.firstFrame;
firstFrame_openSim_restults = round(tempData.time(1)/(1/s.frequency));

if firstFrame_openSim_restults == 0
    firstFrame_openSim_restults_c3d_frame = firstFrame_openSim_restults + s.firstFrame -1;
else
    firstFrame_openSim_restults_c3d_frame = firstFrame_openSim_restults + s.firstFrame -1 - difference_btw_first_frame_and_first_event;
end


plot(frame_first_c3dEvent_in_openSim_frame,'or')

if firstFrame_openSim_restults_c3d_frame > frame_first_c3dEvent_in_openSim_frame
    error(['initial opensim time after first event check file: ' stoFile])
end

preframes = s.preframes + s.firstFrame - firstFrame_openSim_restults_c3d_frame;                    % adjust to match the .sto / .mot file
frame_first_c3dEvent_in_openSim_frame = frame_first_c3dEvent_in_openSim_frame - floor(preframes); % this is a fix for the SO errors --> simulation started a few frames earlier to avoid activation limit


