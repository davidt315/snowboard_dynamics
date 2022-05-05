%
function snowboard
clear all
close all


g = 9.81;                   % gravitational acceleration in m/s^2
m = 50;                     % mass in kg
r = 20;                     % ramp radius in m
COM_board_dist = 0.55*1.754; % distance between snowboard/ground and person's COM. this number is very not fact checked
ramp_theta_cutoff = 15;     % degrees where ramp ends
r_COM = r - COM_board_dist;

% initial values phase 1
z0_1=[0, 0];

% Time and ode
t_span = [0, 15]; 

reltol = 1.0e-6;
options1= odeset('RelTol', reltol,'Events',@event_stop_1);

[T1,Z1] = ode45(@eom1, t_span, z0_1, options1);

T1_end = T1(end);

%Stage 2 setup - defining IC's etc
x0 = r_COM*(1 + sind(ramp_theta_cutoff));%considering left edge of ramp where x = 0
y0 = r_COM*(1 - cosd(ramp_theta_cutoff));%considering lowest point of ramp where y=0
v0_mag = r_COM*(Z1(length(Z1),2));
vx0 = v0_mag*cos(Z1(length(Z1),1) - pi*0.5);
vy0 = v0_mag*sin(Z1(length(Z1),1) - pi*0.5);
z0_2 = [x0,y0,vx0,vy0];
t_span2 = [T1_end, 30];
options2 = odeset('RelTol',reltol,'Events',@event_stop_2);
[T2,Z2] = ode45(@eom2,t_span2,z0_2,options2);

T2 = T2+T1_end;

full_time = [T1;T2];
ramp_angle = Z1(:,1);
ramp_cartesian = [r_COM*(1 - cos(ramp_angle)),r_COM*(1 - sin(ramp_angle))];
air_cartesian = Z2(:,1:2);
cartesian = [ramp_cartesian;air_cartesian];
%{
figure(1)
plot(ramp_cartesian(:,1),ramp_cartesian(:,2));
figure(2)
plot(air_cartesian(:,1),air_cartesian(:,2));
%}
plot(cartesian(:,1),cartesian(:,2));

%options2 = odeset('RelTol',reltol,'Events',event_stop_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = eom1 (T,Z)
    %Stage 1 (Snowboarder on ramp EOM)
% DAE form
% eom + constraint eqn (written wrt state variables) 
%z1 = theta, z2 = thetadot
F_d = 0;%adding a drag force = 0 in case we want to implement this later

dz1dt = Z(2);
dz2dt = g*cos(Z(1))/r - F_d/m;
dzdt = [dz1dt;dz2dt];

%Commenting this version out since it looks like double pendulum stuff
%{
 z1 = theta1, z2 = theta2, z3 = theta1dot, z4 = theta2dot, z5 = T1, z6 =T2

Mat = [1 0 0 0 0 0;
       0 1 0 0 0 0;
       0 0 L1 0 0 ((-sin(Z(2)-Z(1)))/m1);
       0 0 0 (L2*cos(Z(2)-Z(1))) 0 (sin(Z(2)-Z(1))*((1/m1)+(1/m2)));
       0 0 0 0 (-1/m1) (cos(Z(2)-Z(1))/m1);
       0 0 0 (-L2*sin(Z(2)-Z(1))) 0 (cos(Z(2)-Z(1))/m2)];
vec = [Z(3);
       Z(4);
       -g*sin(Z(1));
       (Z(4)^2)*L2*sin(Z(2)-Z(1));
       -(Z(3)^2)*L1 - g*cos(Z(1));
       (Z(3)^2)*L1 + (Z(4)^2)*L2*cos(Z(2)-Z(1)) + g*cos(Z(1))];


dzdt = Mat\vec;
%}
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = eom2(T,Z)
    %Stage 2 - Snowboarder in air's COM EOM
    %z1 = x, z2 = y, z3 = xdot, z4 = ydot
    dzdt = [Z(3);Z(4);0; -g];
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventvalue,stopthecalc,eventdirection] = event_stop_1(T,Z)
    eventvalue      =  (Z(1)-((ramp_theta_cutoff + 90)*2*pi/360));%event at end of ramp defined by ramp angle cutoff
    stopthecalc     =  1;       %  Stop if event occurs
    eventdirection  = 1;       %  Detect only events with dydt>0
end

function [eventval,stop,dir] = event_stop_2(T,Z)
    eventval = Z(2); %event at y = 0
    stop = 1;%stop if event occurs setting = true
    dir = -1;%event only occurs if dydt < 0

end

end