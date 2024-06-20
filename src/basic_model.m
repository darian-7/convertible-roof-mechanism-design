%% ######## Car Roof Model ########
%                09/03/2022 MENG20007
%    Viral Shah, Darian Irani & Andhika Nasution
% ####################################
%% Define System Parameters

clc
clear all
close all

% #### EXTENSIONS ####
% could use interp1 for radius function?

% Mechanism Parameters
g = 9.81; % m/s^2
m = 14.3; % kg [placeholder]
t = 20; % s - time of cycle [placeholder]
c = 1181.102; % damping coefficient
Cd = 2.05; % drag coefficient
rho = 1.225; % density of air
w = 1.204;


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
    
    disp("You have chosen the APM motor. ")
    disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl + "rpm")
    
elseif motorchoice == "1"
    
    nsachoice = input('This category has three motors. Please enter 2, 3 or 4 to select a specific motor ', 's');
    nsachoice = lower(nsachoice);
    
    while nsachoice ~= "2" && nsachoice ~= "3" && nsachoice ~= "4"
    
        motorchoice = input('Please check your input.\n NSA-I [2, 3 or 4]', 's');
        motorchoice = lower(motorchoice);
    
    end
    
    if nsachoice == "2"

        T_s = 0.43;
        w_nl = 2700*2*pi/60;

        disp("You have chosen NSA-I motor one. ")
        disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl + "rad/s")

    elseif nsachoice == "3"

        T_s = 0.48;
        w_nl = 8000*2*pi/60;

        disp("You have chosen NSA-I motor two. ")
        disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl +"rpm")

    elseif nsachoice == "4"

        T_s = 0.69;
        w_nl = 3650*2*pi/60;

        disp("You have chosen NSA-I motor three. ")
        disp("Your selected motor has a stall torque of " + T_s +"Nm " + "and a no load speed of " + w_nl + "rpm")

    end
    
end


    
k = T_s/w_nl; % constant - gradient of line

% Calculate Holding Ratio
r_end = motionpath.r(1)*0.001;
T_gear = abs(m*g*r_end*cos(atan(motionpath.y(1)/motionpath.x(1))));
H_r = T_gear/T_s;
str_H_r = "The minimum gear ratio for this motor is: ";
str_H_r1 = "Please enter 3 gear ratios above this value.";
disp(str_H_r + H_r)
disp(str_H_r1)

% Range of Gear ratios to try
G_r = []; % Work out holding ratio!
G_r(1) = input('What is the first gear ratio? '); % 131.5 or 131 - but the beginning is TOOOOO slow. Have to somehow make it slow down after gravity becomes an accelerating torque.
G_r(2) = input('What is the second gear ratio? ');
G_r(3) = input('What is the third gear ratio? ');
% Create array to store solutions for each tested gear ratio
thet_gr = cell(size(G_r));
time_gr = cell(size(G_r));

%% Define System + ODE


for a0 = 1:length(G_r)
    
    motor = @(t,z) -(k*z(2)*G_r(a0))+m_d*T_s;

    gear = @(t,z) (motor(t,z)*G_r(a0))/(m*radius(t,z,motionpath)^2);

    damping = @(t,z) -(c*z(2))/(m*radius(t,z,motionpath)^2);

    gravity = @(t,z) -g*cos(z(1)-pi/2)/radius(t,z,motionpath); 

    dz = @(t,z)[z(2);
                    gravity(t,z) + gear(t,z) + damping(t,z)];

    IC = [thet_start, 0];
    opts = odeset('RelTol',1e-6,'events', @event_function1);
    [time, theta] = ode45(dz,[0,t],IC,opts);
    thet_gr{a0} = theta;
    time_gr{a0} = time;

    hold on
    set(gcf,'Position',[404 323 1047 421])
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
legend(num2str(G_r(1)), num2str(G_r(2)), num2str(G_r(3)), 'Finishing Angle')


%% ######## Power & Efficiency Calculations ########

% Gear Ratio 1

thet_p1 = thet_gr{1};
time_p1 = time_gr{1};
power1 = zeros(length(thet_gr{1}),1);
torque1 = zeros(length(thet_gr{1}),1);
speedrpm1 = zeros(length(thet_gr{1}),1);

for i1 = 1:length(power1)  
    
    T1 = abs((((-k*thet_p1(i1,2)*G_r(1))+m_d*T_s)));
    P1 = abs(T1*thet_p1(i1,2));
    SRPM1 = abs(thet_p1(i1,2)*G_r(1));
    speedrpm1(i1,1) = SRPM1;
    power1(i1,1) = P1;
    torque1(i1,1) = T1;
    
end



% Gear Ratio 2

thet_p2 = thet_gr{2};
time_p2 = time_gr{2};
power2 = zeros(length(thet_gr{2}),1);
torque2 = zeros(length(thet_gr{2}),1);
speedrpm2 = zeros(length(thet_gr{2}),1);


for i2 = 1:length(power2)  
    
    T2 = abs((((-k*thet_p2(i2,2)*G_r(2))+m_d*T_s)));
    P2 = abs(T2*thet_p2(i2,2));
    SRPM2 = abs(thet_p2(i2,2)*G_r(2));
    speedrpm2(i2,1) = SRPM2;
    power2(i2,1) = P2;
    torque2(i2,1) = T2;
    
end


thet_p3 = thet_gr{3};
time_p3 = time_gr{3};
power3 = zeros(length(thet_gr{3}),1);
torque3 = zeros(length(thet_gr{3}),1);
speedrpm3 = zeros(length(thet_gr{3}),1);

for i3 = 1:length(power3)  
    
    T3 = abs((((-k*thet_p3(i3,2)*G_r(3))+m_d*T_s)));
    P3 = abs(T3*thet_p3(i3,2));
    SRPM3 = abs(thet_p3(i3,2)*G_r(3));
    speedrpm3(i3,1) = SRPM3;
    power3(i3,1) = P3;
    torque3(i3,1) = T3;
    
end

hold off
figure
set(gcf,'Position',[404 323 1047 421])
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

% ######## Torque Speed Relationship ########

% trapz(time_p1,power1) % Energy used in cycle

figure
set(gcf,'Position',[404 323 1047 421])
plot(torque1,speedrpm1)

title('Torque-Speed Curve for Selected Motor: Gear Ratio 1')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Torque [Nm]')
ylabel('Speed [rad/s]')
legend(num2str(G_r(1)))

figure
set(gcf,'Position',[404 323 1047 421])
plot(torque2,speedrpm2)

title('Torque-Speed Curve for Selected Motor: Gear Ratio 2')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Torque [Nm]')
ylabel('Speed [rad/s]')
legend(num2str(G_r(2)))

figure
set(gcf,'Position',[404 323 1047 421])
plot(torque3,speedrpm3)

title('Torque-Speed Curve for Selected Motor: Gear Ratio 3')
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

figure
hold on
plot( thet_p1(:,2), power1)
plot(thet_p2(:,2), power2)
plot(thet_p3(:,2), power3)
hold off

title('Power vs Angular Speed')
lgd = legend;
lgd.Title.String = 'Gear Ratio';
xlabel('Speed [rad/s]')
ylabel('Power [W]')
legend(num2str(G_r(1)),num2str(G_r(2)),num2str(G_r(3)))
