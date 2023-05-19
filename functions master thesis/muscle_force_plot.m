function muscle_force_plot(muscle_force_str, time_str, name, hs1, to, leg)

% plot muscle force
figure('Name',name)
% plot(time_str, muscle_force_stp)
plot(muscle_force_str)
hold on
grid on
% plot(time_stp(1),muscle_force_stp(1),'*k') % first heel strike
% plot(time_stp(to-hs1),muscle_force_stp(to-hs1),'*r') % toe off
% plot(time_stp(end),muscle_force_stp(end),'*k') % second heel strike
plotVert(to); % toe off
% xlabel('gait cycle [%]')
% ylabel('force [BW]')
title(['muscle force of the ' name ' (one step, ' leg ' leg)'])
legend('muscle force', 'toe off')
end

