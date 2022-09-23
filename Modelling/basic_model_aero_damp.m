% ######## Car Roof Model ########
%                09/03/2022 MENG20007
%    Viral Shah, Darian Irani & Andhika Nasution
% ####################################

%% ###### Define System Parameters ######

clc
clear all
close all

% Mechanism Parameters
g = 9.81; % m/s^2
snow = 11.26;
m = 23.5; % kg [placeholder]
t = 20; % s - time of cycle [placeholder]
Cd = 2.05; % drag coefficient
rho = 1.225; % density of air
% w = 1.204;
windmph = 20; % target from PDS
wind = windmph/2.237;
r_j = 0.1031; % distance from damper connector to motor [m]
R = 0.370;

% ## DAMPING ## 

% Linear Bi-directional
l1 = 5061.867;
l2 = 7086.614; 
l3 = 11811.024;
l4 = 17716.535;

% Rotational Damper Coefficient Nm/rad/s

crd = input('What is the damping coefficient in Nm/rad/s? ');


%% ###### Deploy or Retract ######

dor = input('Deploy or Retract? [0 or 1] ', 's');
dor = lower(dor);
load('com.mat');

if dor == "0"
    
    thet_start = pi/2 + 0.4186630; % rad - start angle [0.4186632]
    thet_end = pi/2 + 2.6213734; % rad - finish angle [2.6203882]
    m_d = 1;
    
elseif dor == "1"
    
    thet_end = pi/2 + 0.4186630; % rad - start angle [0.4186632]
    thet_start = pi/2 + 2.6213734; % rad - finish angle [2.6203882]
    m_d = -1;
    
end
    

while dor ~= "1" && dor ~= "0"
    
    dor = input('Please check your input.\nDeploy or Retract? [0 or 1] ', 's');
    dor = lower(dor);
    
    if dor == "0"
    
        thet_start = pi/2 + 0.4186630; % rad - start angle [0.4186632]
        thet_end = pi/2 + 2.6213734; % rad - finish angle [2.6203882]
        m_d = 1;

    elseif dor == "1"

        thet_end = pi/2 + 0.4186630; % rad - start angle [0.4186632]
        thet_start = pi/2 + 2.6213734; % rad - finish angle [2.6203882]
        m_d = -1;

    end
    
end

%% ####### Selection of Motors #######

% ######## Motor Choices ########

motorchoice = input('What motor to choose? [type 0 for APM / type 1 for NSA-I]: ', 's');
motorchoice = lower(motorchoice);

while motorchoice ~= "1" && motorchoice ~= "0"
    
    motorchoice = input('Please check your input.\nAPM or NSA-1? [0 or 1] ', 's');
    motorchoice = lower(motorchoice);
    
end

if motorchoice == "0"
    
    T_s = 0.19;
    w_nl = 5000*2*pi/60;
    mt = 43.63;
    c = 0.66;
    
    disp("You have chosen the APM motor. ")
    disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl + "rad/s")
    
end
    
if motorchoice == "1"
    
    nsachoice = input('NSA-I 1, 2 or 3? ', 's');
    nsachoice = lower(nsachoice);
    
    while nsachoice ~= "1" && nsachoice ~= "2" && nsachoice ~= "3"

    nsachoice = input('Please check your input.\n NSA-I [1, 2 or 3]', 's');
    nsachoice = lower(nsachoice);

    end
    
    if nsachoice == "1"

        T_s = 0.43;
        w_nl = 2700*2*pi/60;
        mt = 30.16;
        c = 1.09;

        disp("You have chosen NSA-I motor one. ")
        disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl + "rad/s")

    elseif nsachoice == "2"

        T_s = 0.48;
        w_nl = 8000*2*pi/60;
        mt = 77.4;
        c = 4.75;

        disp("You have chosen NSA-I motor two. ")
        disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl +"rad/s")

    elseif nsachoice == "3"

        T_s = 0.69;
        w_nl = 3650*2*pi/60;
        mt = 43.48;
        c = 1.43;

        disp("You have chosen NSA-I motor three. ")
        disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl + "rad/s")

    end
    
end
    
   
k = T_s/w_nl; % constant - gradient of line

% Calculate Holding Ratio
r_end = motionpath.r(1)*0.001;
T_gear = abs(m*g*r_end*cos(motionpath.thet(1)+pi/2));
H_r = T_gear/T_s;
str_H_r = "The minimum gear ratio for this motor is: ";
str_H_r1 = "Please enter 3 gear ratios above this value.";
disp(str_H_r + H_r)
disp(str_H_r1)

% Range of Gear ratios to try
G_r = []; % Work out holding ratio!
G_r(1) = input('What is the first gear ratio? '); 
G_r(2) = input('What is the second gear ratio? ');
G_r(3) = input('What is the third gear ratio? ');
% Create array to store solutions for each tested gear ratio
thet_gr = cell(size(G_r));
time_gr = cell(size(G_r));

%% ######## Define System + ODE ########


for a0 = 1:length(G_r)
    
    ad = @(t,z) (0.5*rho*Cd*aero(t,z,motionpath)*(wind^2)*sin(z(1)-pi/2))/(m*radius(t,z,motionpath));
    
    motor = @(t,z) -(k*z(2)*G_r(a0))+m_d*T_s;

    gear = @(t,z) (motor(t,z)*G_r(a0))/(m*radius(t,z,motionpath)^2);

%     gear = @(t,z) (motor(t,z)*G_r(a0))/(m*R^2);

    damper = @(t,z) -damping(t,z,motionpath,m_d)*(z(2))*r_j^2/(m*radius(t,z,motionpath)^2);
    
    % LINEAR DAMPER
    
    Ldamper = @(t,z) -l4*(z(2))*r_j^2/(m*radius(t,z,motionpath)^2);
    
    % DASHPOT
    
    Rdamper = @(t,z) -crd*(z(2))/(m*radius(t,z,motionpath)^2);
    
    gravity = @(t,z) -g*cos(z(1)-pi/2)/radius(t,z,motionpath); 
    
%     gravity = @(t,z) -g*cos(z(1)-pi/2)/R; 
    
    dz = @(t,z)[z(2);
                    gravity(t,z) + gear(t,z) + Rdamper(t,z) - ad(t,z)];
             
%     dz = @(t,z)[z(2);
%                     gravity(t,z) + gear(t,z)];

    IC = [thet_start, 0];
    opts = odeset('RelTol',1e-7,'events', @event_function1);
    [time, theta] = ode45(dz,[0,t],IC,opts);
    thet_gr{a0} = theta;
    time_gr{a0} = time;

    hold on
%     set(gcf,'Position',[404 323 1047 421])
    plot(time,theta(:,1))
%     animate_pendulum(time,theta(:,1));
    disp(time(end))

end

if dor == "0"
    
    title('Angle vs Time for Deploying Mechanism')
    
elseif dor == "1"
    
    title('Angle vs Time for Retracting Mechanism')
    
end

lgd = legend; % Position: 0.1369,0.6956,0.128,0.2126
lgd.Title.String = 'Gear Ratio';
xlabel('Time [s]')
ylabel('\theta [rad]')
yline(thet_end)
xline(15)
xline(20)
legend(num2str(G_r(1)), num2str(G_r(2)), num2str(G_r(3)), 'Finishing Angle')
% legend(num2str(G_r(1)), 'Finishing Angle')
% saveas(gcf, 'Angle vs Time.png')
%% ######## Output Graphs ########

% Gear Ratio 1

thet_p1 = thet_gr{1};
time_p1 = time_gr{1};
damping1 = zeros(length(thet_gr{1}),1);
power1 = zeros(length(thet_gr{1}),1);
torque1 = zeros(length(thet_gr{1}),1);
motorspeed1 = zeros(length(thet_gr{1}),1);
current1 = zeros(length(thet_gr{1}),1);
powerin1 = zeros(length(thet_gr{1}),1);
eff1 = zeros(length(thet_gr{1}),1);

for i1 = 1:length(power1)  
    
    T1 = abs((((-k*thet_p1(i1,2)*G_r(1))+m_d*T_s)));
    P1 = abs(T1*thet_p1(i1,2)*G_r(1)); % output power
    D1 = crd*(thet_p1(i1,2));
    I1 = mt*T1 + c;
    p1 = 12*I1;
    MS1 = abs(thet_p1(i1,2)*G_r(1));
    e1 = 100*P1/p1;
    eff1(i1,1) = e1;
    motorspeed1(i1,1) = MS1;
    power1(i1,1) = P1;
    torque1(i1,1) = T1;
    damping1(i1,1) = D1;
    current1(i1,1) = I1;
    powerin1(i1,1) = p1;
    
end



% Gear Ratio 2

thet_p2 = thet_gr{2};
time_p2 = time_gr{2};
damping2 = zeros(length(thet_gr{2}),1);
power2 = zeros(length(thet_gr{2}),1);
torque2 = zeros(length(thet_gr{2}),1);
motorspeed2 = zeros(length(thet_gr{2}),1);
current2 = zeros(length(thet_gr{2}),1);
powerin2 = zeros(length(thet_gr{2}),1);
eff2 = zeros(length(thet_gr{2}),1);


for i2 = 1:length(power2)  
    
    T2 = abs((((-k*thet_p2(i2,2)*G_r(2))+m_d*T_s)));
    P2 = abs(T2*thet_p2(i2,2)*G_r(2));
    D2 = crd*(thet_p2(i2,2));
    I2 = mt*T2 + c;
    p2 = 12*I2;
    MS2 = abs(thet_p2(i2,2)*G_r(2));
    e2 = 100*P2/p2;
    eff2(i2,1) = e2;
    motorspeed2(i2,1) = MS2;
    power2(i2,1) = P2;
    torque2(i2,1) = T2;
    damping2(i2,1) = D2;
    current2(i2,1) = I2;
    powerin2(i2,1) = p2;
    
end

% Gear Ratio 1

thet_p3 = thet_gr{3};
time_p3 = time_gr{3};
damping3 = zeros(length(thet_gr{3}),1);
power3 = zeros(length(thet_gr{3}),1);
torque3 = zeros(length(thet_gr{3}),1);
motorspeed3 = zeros(length(thet_gr{3}),1);
current3 = zeros(length(thet_gr{3}),1);
powerin3 = zeros(length(thet_gr{3}),1);
eff3 = zeros(length(thet_gr{3}),1);

for i3 = 1:length(power3)  
    
    T3 = abs((((-k*thet_p3(i3,2)*G_r(3))+m_d*T_s)));
    P3 = abs(T3*thet_p3(i3,2)*G_r(3));
    D3 = crd*(thet_p3(i3,2));
    I3 = mt*T3 + c;
    p3 = 12*I3;
    MS3 = abs(thet_p3(i3,2)*G_r(3));
    e3 = 100*P3/p3;
    eff3(i3,1) = e3;
    motorspeed3(i3,1) = MS3;
    power3(i3,1) = P3;
    torque3(i3,1) = T3;
    damping3(i3,1) = D3;
    current3(i3,1) = I3;
    powerin3(i3,1) = p3;
    
end

figure
% set(gcf,'Position',[404 323 1047 421])
hold on
plot(time_p1,power1)
plot(time_p2,power2)
plot(time_p3,power3)
hold off

title('Mechanical Power of Mechanism vs Time')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Time [s]')
ylabel('Power [W]')
legend(num2str(G_r(1)), num2str(G_r(2)), num2str(G_r(3)))
% legend(num2str(G_r(1)))
% saveas(gcf, 'Mechanical Power vs Time.png')

% #### Energy Usage ####

cn = newline;
E1 = trapz(time_p1,power1); % Energy used in cycle
E2 = trapz(time_p2,power2);
E3 = trapz(time_p3,power3);

E_S = "The energy used with each respective gear ratio is as follows: ";

disp(E_S + cn + E1+"J" + cn + E2+"J" +  cn + E3+"J")
% disp(E_S + E1+"J")

EE1 = trapz(time_p1,powerin1); % Energy used in cycle
EE2 = trapz(time_p2,powerin2);
EE3 = trapz(time_p3,powerin3);

EE_S = "The energy supplied with each respective gear ratio is as follows: ";

disp(cn + EE_S + cn + EE1+"J" + cn + EE2+"J" +  cn + EE3+"J")
% disp(E_S + E1+"J")

EFF1 = E1*100/EE1;
EFF2 = E2*100/EE2;
EFF3 = E3*100/EE3;

disp("Efficiencys are:" + cn + EFF1+"%" + cn + EFF2+"%" +  cn + EFF3+"%")

figure
hold on
% set(gcf,'Position',[404 323 1047 421])
subplot(1,3,1)
plot(torque1,motorspeed1)

title('Torque vs Angular Speed @ Motor')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Torque [Nm]')
ylabel('Speed [rad/s]')
legend(num2str(G_r(1)))

subplot(1,3,2)
plot(torque2,motorspeed2)

title('Torque vs Angular Speed @ Motor')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Torque [Nm]')
ylabel('Speed [rad/s]')
legend(num2str(G_r(2)))

subplot(1,3,3)
plot(torque3,motorspeed3)
hold off

title('Torque vs Angular Speed @ Motor')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Torque [Nm]')
ylabel('Speed [rad/s]')
legend(num2str(G_r(3)))

figure
hold on
plot(time_p1, abs(thet_p1(:,2)))
plot(time_p2, abs(thet_p2(:,2)))
plot(time_p3, abs(thet_p3(:,2)))
hold off

title('Angular Speed vs Time')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Time [s]')
ylabel('Speed [rad/s]')
legend(num2str(G_r(1)),num2str(G_r(2)),num2str(G_r(3)))
% legend(num2str(G_r(1)))
% saveas(gcf, 'Angular Speed vs Time.png')
figure
hold on
plot(thet_p1(:,2), power1)
plot(thet_p2(:,2), power2)
plot(thet_p3(:,2), power3)
hold off

title('Power vs Angular Speed')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Speed [rad/s]')
ylabel('Power [W]')
legend(num2str(G_r(1)),num2str(G_r(2)),num2str(G_r(3)))
% legend(num2str(G_r(1)))
% saveas(gcf, 'Power vs Angular Speed.png')

figure
% set(gcf,'Position',[404 323 1047 421])
hold on
plot(time_p1,abs(damping1))
plot(time_p2,abs(damping2))
plot(time_p3,abs(damping3))
hold off

title('Damping Torque vs Time')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Time [s]]')
ylabel('Damping Torque [Nm]')
legend(num2str(G_r(1)),num2str(G_r(2)),num2str(G_r(3)))

% figure
% hold on
% plot(torque1,thet_p1(:,2))
% plot(torque2,thet_p2(:,2))
% plot(torque3,thet_p3(:,2))
% hold off

figure
hold on
plot(time_p1,eff1)
plot(time_p2,eff2)
plot(time_p3,eff3)
hold off

title('Motor Efficiency vs Time')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Time [s]')
ylabel('Motor Efficiency [%]')
legend(num2str(G_r(1)),num2str(G_r(2)),num2str(G_r(3)))

disp('The max power for each gear ratio is: ')
disp(max(power1))
disp(max(power2))
disp(max(power3))

