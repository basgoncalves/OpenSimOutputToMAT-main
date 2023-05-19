function [moments_stp] = moments_plot(moments_stp,time_stp, name ,hs1, to, leg)

% plot moments
figure('Name',name)
plot(time_stp, moments_stp)
hold on
grid on
%plot(time_stp(1),moments_stp(1),'*k') % first heel strike
plot(time_stp(to-hs1),moments_stp(to-hs1),'*r') % toe off
%plot(time_stp(end),moments_stp(end),'*k') % second heel strike
xlabel('gait cycle [%]')
ylabel('moment [BW]')
title([name ' moment (one step, ' leg ' leg)'])
legend('moment','toe off')
end

