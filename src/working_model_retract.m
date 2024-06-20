%% Car Roof Model

% ####################################

% 09/03/2022 MENG20007
% Viral Shah, Darian Irani & Andhika Nasution

% Possible Extensions: drag, damping, variable mass/radius
% VARIABLE RADIUS function of z(1)...

% ####################################
%% Define System Parameters

clc
clear all

dor = input('Deploy or Retract? [1 or 0] ', 's');
dor = lower(dor);
load('com.mat');

while dor ~= "1" && dor ~= "0"
    
    dor = input('Please check your spelling.\nDeploy or Retract? [1 or 0] ', 's');
    dor = lower(dor);
    
end

% Mechanism Parameters
g = 9.81; % m/s^2
m = 20; % kg [placeholder]
t = 15; % s - time of cycle [placeholder]
c = [75]; % damping coefficient

thet_in = pi/2 + 0.4186632; % rad - start angle [0.4186632]
thet_out = pi/2 + 2.6203882; % rad - finish angle [2.6203882]

% Motor Constants
T_s = 0.69; % Nm - stall torque
% w_m = 4000*2*pi/60; % rad/s - motor speed
w_nl = 3650*2*pi/60; % rad/s - no-load speed
k = T_s/w_nl; % constant - gradient of line

% Range of Gear ratios to try
G_r = [140, 150, 160]; % Work out holding ratio!

% Create array to store solutions for each tested gear ratio
thet_gr = cell(size(G_r));

%% Define System + ODE

if dor == "1"

    for a0 = 1:length(G_r)
    
        for a1 = 1:length(c)

            motor = @(t,z) ((-(k*z(2)*G_r(a0))+T_s)*G_r(a0))/(m*radius(t,z,motionpath)^2);

            damping = @(t,z) -c(a1)*z(2)/(m*radius(t,z,motionpath)^2);

            gravity = @(t,z) -g*cos(z(1)-pi/2)/radius(t,z,motionpath); % beyond 90 degrees, gravity = accelerating torque: is this included already?

            dz = @(t,z)[z(2);
                            gravity(t,z) + motor(t,z) + damping(t,z)];


            IC = [thet_in, 0];
            opts = odeset('RelTol',1e-6,'events', @event_function1);
            [time, theta] = ode45(dz,[0,t],IC,opts);
            thet_gr{a0} = theta;

            plot(time,theta(:,1))
            xlabel('time [s]')
            ylabel('rad [theta]') 
            hold on
            yline(thet_out)
%             animate_pendulum(time,theta(:,1));
            hold on

        end
        
    end

end


% Retract not working - how to make motor go opposite direction? Tried making accel. torque negative but breaks the script.


if dor == "0"

    for a0 = 1:length(G_r)
    
        for a1 = 1:length(c)

            motor = @(t,z) ((-(k*z(2)*G_r(a0))-T_s)*G_r(a0))/(m*radius(t,z,motionpath)^2);

            damping = @(t,z) -c(a1)*z(2)/(m*radius(t,z,motionpath)^2);

            gravity = @(t,z) -g*cos(z(1)-pi/2)/radius(t,z,motionpath); % beyond 90 degrees, gravity = accelerating torque: is this included already?

            dz = @(t,z)[z(2);
                            gravity(t,z) + motor(t,z) + damping(t,z)];


            IC = [thet_out, 0];
            opts = odeset('RelTol',1e-6,'events', @event_function2);
            [time, theta] = ode45(dz,[0,t],IC,opts);
            thet_gr{a0} = theta;

            plot(time,theta(:,1))
            xlabel('time [s]')
            ylabel('rad [theta]') 
            hold on
            yline(thet_in)
%             animate_pendulum(time,theta(:,1));
            hold on

        end
        
    end

end



