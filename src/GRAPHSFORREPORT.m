clear all
close all

%% Technical Analysis First 3 Extensions

thet_endR = pi/2 + 0.4186630;
thet_endD = pi/2 + 2.6213734;

figure

hold on

load('CRthetaR.mat');
load('CRtimeR.mat');

plot(time,theta(:,1), 'k') % black

load('VRthetaR.mat');
load('VRtimeR.mat');

plot(time,theta(:,1),'k--') % dashed

load('AEROthetaR.mat');
load('AEROtimeR.mat');

plot(time,theta(:,1),'k-.') % 

load('CRthetaD.mat');
load('CRtimeD.mat');

plot(time_p1,thet_p1(:,1), 'b') % 

load('VRthetaD.mat');
load('VRtimeD.mat');

plot(time_p1,thet_p1(:,1),'b--') % dashed

load('AEROthetaD.mat');
load('AEROtimeD.mat');

plot(time_p1,thet_p1(:,1),'b-.') % 

hold off

title('Angle vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = ('Retracting [black] & Deploying [blue]');
xlabel('Time [s]','FontSize',12)
ylabel('\theta [rad]','FontSize',12)
yline(thet_endR,'k:')
yline(thet_endD,'b:')
legend('Constant Radius', 'Variable Radius', 'Aerodynamics','Constant Radius', 'Variable Radius', 'Aerodynamics', 'Finishing Angle', 'Finishing Angle','NumColumns',2)
set(gca, 'FontName', 'Times New Roman')

%% Damping Coefficient for same GR

figure
hold on

load('DAMP1time.mat');
load('DAMP1angle.mat');

plot(time_p1,thet_p1(:,1), 'k') % black

load('DAMP2time.mat');
load('DAMP2angle.mat');

plot(time_p1,thet_p1(:,1),'k--') % dashed

load('DAMP3time.mat');
load('DAMP3angle.mat');

plot(time_p1,thet_p1(:,1),'k-.') % blue

load('DAMP600time.mat');
load('DAMP600angle.mat');

plot(time_p1,thet_p1(:,1),'b') % dashed

load('DAMP1000time.mat');
load('DAMP1000angle.mat');

plot(time_p1,thet_p1(:,1),'b--') % blue

hold off

title('Angle vs Time')
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Dampers & Coefficient';
xlabel('Time [s]','FontSize',12)
ylabel('\theta [rad]','FontSize',12)
yline(thet_endR, 'k:')
legend('Extension: 7086', 'Extension: 17716', 'Bi-directional: 17716', 'CRD: 600', 'CRD: 100', 'Finishing Angle')
set(gca, 'FontName', 'Times New Roman')

%% Damping Torque Speed Time Graph

load('DAMPfinaltime');
load('DAMPfinaltorque');
load('DAMPfinalspeed');

figure
yyaxis left
plot(time_p1,abs(damping1./421.45),'k')
ylabel('Damping Torque [Nm]','FontSize',12)
yyaxis right
plot(time_p1,abs(thet_p1(:,2).*421.45),'k--')
ylabel('Angular Speed [rad/s]','FontSize',12)
hold off

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

title('Damping & Angular Speed @ Motor vs Time','FontSize',12)
xlabel('Time [s]','FontSize',12)
legend('Torque', 'Angular Speed')
set(gca, 'FontName', 'Times New Roman')
%% Initial Model Motor

load('initialAPMangle');
load('initialAPMtime');

figure
hold on
plot(time_p1,thet_p1(:,1), 'k') % black

load('initialNSA1angle');
load('initialNSA1time');

plot(time_p1,thet_p1(:,1), 'k--') 

load('initialNSA2angle');
load('initialNSA2time');

plot(time_p1,thet_p1(:,1), 'k:') 

load('initialNSA3angle');
load('initialNSA3time');

plot(time_p1,thet_p1(:,1), 'k-.') 
hold off

title('Angle vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('\theta [rad]','FontSize',12)
yline(thet_endR)
legend('APM', 'NSA-I 1', 'NSA-I 2', 'NSA-I 3')
set(gca, 'FontName', 'Times New Roman')
%% Initial Model Efficiency

load('initialAPMtime');
load('initialAPMeff');

figure
hold on
plot(time_p1,eff1, 'k') % black

load('initialNSA1eff');
load('initialNSA1time');

plot(time_p1,eff1, 'k--') 

load('initialNSA2eff');
load('initialNSA2time');

plot(time_p1,eff1, 'k:') 

load('initialNSA3eff');
load('initialNSA3time');

plot(time_p1,eff1, 'k-.') 
hold off

title('Efficiency vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('Efficiency [%]','FontSize',12)
legend('APM', 'NSA-I 1', 'NSA-I 2', 'NSA-I 3')
set(gca, 'FontName', 'Times New Roman')

%% Final Model Motor

% ANGLE

figure
hold on

load('finalNSA1angle');
load('finalNSA1time');

plot(time_p1,thet_p1(:,1), 'k--') 

load('finalNSA2angle');
load('finalNSA2time');

plot(time_p1,thet_p1(:,1), 'k:') 

load('finalNSA3angle');
load('finalNSA3time');

plot(time_p1,thet_p1(:,1), 'k') 
hold off

title('Angle vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('\theta [rad]','FontSize',12)
yline(thet_endR)
legend('NSA-I 1', 'NSA-I 2', 'NSA-I 3')
set(gca, 'FontName', 'Times New Roman')

% POWER

figure
hold on

load('finalNSA1power');
load('finalNSA1time');

plot(time_p1,power1, 'k--') 

load('finalNSA2power');
load('finalNSA2time');

plot(time_p1,power1, 'k:') 

load('finalNSA3power');
load('finalNSA3time');

plot(time_p1,power1, 'k') 
hold off

title('Power vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('Power [W]','FontSize',12)
legend('NSA-I 1', 'NSA-I 2', 'NSA-I 3')
set(gca, 'FontName', 'Times New Roman')

%% Final Model Efficiency

figure
hold on

load('finalNSA1eff');
load('finalNSA1time');

plot(time_p1,eff1, 'k--') 

load('finalNSA2eff');
load('finalNSA2time');

plot(time_p1,eff1, 'k:') 

load('finalNSA3eff');
load('finalNSA3time');

plot(time_p1,eff1, 'k') 
hold off

title('Efficiency vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('Efficiency [%]','FontSize',12)
legend('NSA-I 1', 'NSA-I 2', 'NSA-I 3')
set(gca, 'FontName', 'Times New Roman')

%% FINAL REVIEW GRAPHS

figure
hold on

load('FINALtimeR');
load('FINALangleR');

plot(time_p1,thet_p1(:,1), 'k') % black

load('FINALtimeD');
load('FINALangleD');

plot(time_p1,thet_p1(:,1), 'b') 

load('FINALsnowtimeR');
load('FINALsnowangleR');

plot(time_p1,thet_p1(:,1), 'k--') % black

load('FINALtimesnowD');
load('FINALanglesnowD');

plot(time_p1,thet_p1(:,1), 'b--') 

hold off

title('Angle vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = ('Retracting [black] & Deploying [blue]');
xlabel('Time [s]','FontSize',12)
ylabel('\theta [rad]','FontSize',12)
yline(thet_endR,'k:')
yline(thet_endD,'b:')
legend('Retraction', 'Deploy','Retraction Snow', 'Deploy Snow','Finishing Angle','Finishing Angle','NumColumns',2)
set(gca, 'FontName', 'Times New Roman')
load('FINALtimeR');
load('FINALangleR');

figure
hold on

load('FINALtimeR');
load('FINALpowerR');

plot(time_p1,power1, 'k') % black

load('FINALtimeD');
load('FINALpowerD');

plot(time_p1,power1, 'k--') 
hold off

title('Power vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('Power [W]','FontSize',12)
legend('Retraction', 'Deploy')
set(gca, 'FontName', 'Times New Roman')

figure
hold on

load('FINALtimeR');
load('FINALeffR');

plot(time_p1,eff1, 'k') % black

load('FINALtimeD');
load('FINALeffD');

plot(time_p1,eff1, 'k--') 
hold off

title('Efficiency vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor';
xlabel('Time [s]','FontSize',12)
ylabel('Efficiency [%]','FontSize',12)
legend('Retraction', 'Deploy')
set(gca, 'FontName', 'Times New Roman')

%% Damping Coefficient Selection 

figure
hold on

load('350NSA1angle');
load('350NSA1time');

plot(time_p1,thet_p1(:,1), 'k:') % black

load('1350NSA1angle');
load('1350NSA1time');

plot(time_p1,thet_p1(:,1), 'k--') % black --

load('1600NSA2angle');
load('1600NSA2time');

plot(time_p1,thet_p1(:,1), 'k-.') % blue

load('2100NSA3angle');
load('2100NSA3time');

plot(time_p1,thet_p1(:,1), 'k') % blue :
hold off

title('Angle vs Time','FontSize',12)
lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Motor/Gear Ratio/Damping Coefficient';
xlabel('Time [s]','FontSize',12)
ylabel('\theta [rad]','FontSize',12)
yline(thet_endR)
legend('NSA-I 1 / 280 / 350', 'NSA-I 1 / 650 / 1350', 'NSA-I 2 / 510 / 1600', 'NSA-I 3 / 450 / 2100')
set(gca, 'FontName', 'Times New Roman')

%% POWER SPEED

% figure
% hold on
% 
% load('420speed');
% load('420power');
% 
% plot(motorspeed1,power1, 'k') % black
% 
% load('435speed');
% load('435power');
% 
% plot(motorspeed2,power2, 'k') % black
% 
% load('450speed');
% load('450power');
% 
% plot(motorspeed3,power3, 'k') % black
% 
% hold off
% 
% title('Power vs Motor Speed','FontSize',12)
% lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
% lgd.Title.String = 'Motor/Gear Ratio/Damping Coefficient';
% xlabel('Time [s]','FontSize',12)
% ylabel('\theta [rad]','FontSize',12)
% yline(thet_endR)
% legend('NSA-I 1 / 280 / 350', 'NSA-I 1 / 650 / 1350', 'NSA-I 2 / 510 / 1600', 'NSA-I 3 / 450 / 2100')
% set(gca, 'FontName', 'Times New Roman')

