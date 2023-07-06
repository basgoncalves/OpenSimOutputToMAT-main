
%% modify figures
% for i=1:4
%     f.Children(i).FontSize = 15;
%     if i>1
%      f.Children(i).XTickLabelRotation = 90;
%     end
% end
% f.Children(3).FontSize = 15;
% f.Children(3).FontSize = 20;

%% create figure GRF
addpath(genpath('C:\Users\Balu\Nextcloud\Documents\MA\Code\MSKmodelling'))
grf = importdata('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD04\C3D\3DGAIT_B_W9\grf.mot');
plot(grf.data(:,3))
xlabel('Time Frames')
ylabel('Vertical Ground Reaction Force [N]')
f=gcf;
f.FontName = 'Arial';
f.FontSize = 10;