% gather participants' characteristics
clc; %clear all; close all

part_char = struct;

% pre surgery
part_char.session1.left.NSA = [128.9, 134.1, 140.3, 137.0, 139.3]; 
part_char.session1.left.AVA = [42.2, 27.2, 8.0, 36.3, 12.7];
part_char.session1.left.TT = [0.1, 6.3, 49.4, 37.7, 45.8];
part_char.session1.right.NSA = [135.2, 129.0, 139.7, 130.6, 145.5];
part_char.session1.right.AVA = [52.3, 22.2, 14.6, 42.1, 18.0];
part_char.session1.right.TT = [3.4, 8.7, 34.1, 50.2, 47.3];
part_char.session1.mass = [32.9, 39.3, 55.9, 40.9, 66.3];


% post surgery
part_char.session2.left.NSA = part_char.session1.left.NSA;
part_char.session2.left.AVA = part_char.session1.left.AVA;
part_char.session2.left.TT = [10.1, 21.3, 19.4, 17.7, 27.8];
part_char.session2.right.NSA = part_char.session1.right.NSA;
part_char.session2.right.AVA = part_char.session1.right.AVA;
part_char.session2.right.TT = [18.4, 23.7, 14.1, 25.2, 25.3];
part_char.session2.mass = [62.5, 39.1, 55.7, 49.7, 72.4];

% TDs
part_char.session3.left.NSA = [121.1, 120.1, 126.2, 125.8];
part_char.session3.left.AVA = [19.8, 27.3, 20.2, 26.1];
part_char.session3.left.TT = [23.1, 25.8, 8.7, 3.4];
part_char.session3.right.NSA = [132.4, 126, 124.6, 124.8];
part_char.session3.right.AVA = [24.2, 13.3, 17.7, 28.7];
part_char.session3.right.TT = [13.7, 26.9, 22.2, 4.9];
part_char.session3.mass = [38.9, 54.8, 29.6, 20.4];

% determine walking velocity
f = fields(data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final);
pel_vel = [];
for i=1:length(f)
    pel_x = data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final.(f{i}).pelvis_tx;
    pel_time = data.IK.FINAL_PERSONALISEDTORSIONS_scaled_final.(f{i}).time;
    pel_dist = pel_x(end)-pel_x(1);
    pel_t = pel_time(end) - pel_time(1);
    pel_vel(i) = pel_dist/pel_t;
end



% save data
save C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\participants_characteristics.mat part_char -v7.3