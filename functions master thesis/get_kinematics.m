function [kin] = get_kinematics(mot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
kin.time = mot.data(:,1);
kin.hip_fl_r = mot.data(:,8);
kin.hip_ad_r = mot.data(:,9);
kin.hip_rot_r= mot.data(:,10);
kin.knee_angle_r = mot.data(:,11);
kin.ankle_angle_r = mot.data(:,12);
kin.hip_fl_l = mot.data(:,15);
kin.hip_ad_l = mot.data(:,16);
kin.hip_rot_l = mot.data(:,17);
kin.knee_angle_l = mot.data(:,18);
kin.ankle_angle_l = mot.data(:,19);
end