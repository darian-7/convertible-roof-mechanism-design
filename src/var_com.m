%% Variable Centre of Mass
clear all

% Graphs
com = readtable('motion_path.csv');
subplot(1,2,1)
set(gcf,'Position',[404 323 1047 421])
plot(com.x,com.y)
title('Position of Centre of Mass for Mechanism')
xlabel('X-coordinate [mm]')
ylabel('Y-coordinate [mm]')
hold on
subplot(1,2,2)
plot(com.thet,com.r)
title('Varying Radius against Angle of Mechanism')
xlabel('Theta [rad]')
ylabel('Radius [mm]')
