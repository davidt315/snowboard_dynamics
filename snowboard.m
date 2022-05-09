%
function snowboard
clear all
close all


g = 9.81;                   % gravitational acceleration in m/s^2
m = 50;                     % mass in kg
r = 20;                     % ramp radius in m
ramp_theta_cutoff = 50;     % where ramp ends, degrees


%Model the snowboarder as a cylinder rotating around an axis through the
%"long" side, not the circular faces. The snowboarder has an MOI of 1/12
%ml^2, where l = the distance from their COM to the board. We change the
%length of the cylinder at a constant rate.
COM_board_dist = 0.55*1.754; % distance btwn snowboard/ground and person's COM. this number is very not fact checked
COM_board_end = 0.3048;     % distance btwn snowboard/ground and person's COM, crouched
MOI_change_time = 3;      % time taken for moment of inertia to change
Ldot = (COM_board_end - COM_board_dist)/MOI_change_time;
r_COM = r - COM_board_dist;

% initial values phase 1
z0_1=[0, 0];

% Time and ode
t_span = [0, 15]; 

reltol = 1.0e-9;
options1= odeset('RelTol', reltol,'Events',@event_stop_1);

[T1,Z1] = ode45(@eom1, t_span, z0_1, options1);

T1_end = T1(end);

%Stage 2 setup - defining IC's etc
x0 = r + r_COM*sind(ramp_theta_cutoff);%considering left edge of ramp where x = 0
y0 = r - r_COM*cosd(ramp_theta_cutoff);%considering lowest point of ramp where y=0
v0_mag = r_COM*(Z1(length(Z1),2));
vx0 = v0_mag*cos(Z1(length(Z1),1) - pi*0.5);
vy0 = v0_mag*sin(Z1(length(Z1),1) - pi*0.5);
%angle of snowboarder's body - COM relative to edge
t_b0 = Z1(length(Z1),1);
t_bdot0 = Z1(length(Z1),2);
z0_2 = [x0,y0,vx0,vy0,t_b0,t_bdot0,COM_board_dist];
t_span2 = [T1_end, 30];
options2 = odeset('RelTol',reltol,'Events',@event_stop_3);

[T2,Z2] = ode45(@eom2,t_span2,z0_2,options2);
T2 = T2+T1_end;
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%This was an attempt to use a third ODE to
model the constant MOI phase of the motion. I realized I couldn't use the
events on it properly, so it's commented out

%Stage 3 setup - defining ICs
end_index = length(Z2);
T2_end = T2(end);
z0_3 = Z2(end_index,1:6);
t_span3 = [T2_end,40];
options3 = odeset('RelTol',reltol,'Events',@event_stop_3);
[T3,Z3] = ode45(@eom3,t_span3,z0_3,options3);
%}
full_time = [T1;T2];
ramp_angle = Z1(:,1);
%Convert stage 1 to Cartesian Coordinates
COMramp_cart = [r - r_COM*cos(ramp_angle),r - r_COM*sin(ramp_angle)];
Edge_ramp_cart = [r*(1 - cos(ramp_angle)),r*(1 - sin(ramp_angle))];

%Get stage 2 position in Cartesian Coordinates
COMair_cart2 = Z2(:,1:2);
body_angle_2 = Z2(:,5);
Edge_air_cart2 = COMair_cart2 - (Z2(:,7)).*[cos(body_angle_2),sin(body_angle_2)];


%{
%%%%%%%%%%%%%%%% Removed stage 3 stuff
COMair_cart3 = Z3(:,1:2);
body_angle_3 = Z3(:,5);
Edge_air_cart3 = COMair_cart3 - COM_board_end.*[cos(body_angle_3),sin(body_angle_3)];
Edge_air_cart = [Edge_air_cart2;Edge_air_cart3];
%}
%Combine coordinates
COMcart = [COMramp_cart;COMair_cart2];
Edge_cart = [Edge_ramp_cart;Edge_air_cart2];

%{
figure(1)
%%%%%%%%%%%%%%% Plotting stage 1 and 2 separately for debugging, no longer
necessary for final 
plot(COMramp_cart(:,1),COMramp_cart(:,2));
figure(2)
plot(COMair_cart2(:,1),COMair_cart2(:,2));
%}
hold on
plot(COMcart(:,1),COMcart(:,2));
plot(Edge_cart(:,1),Edge_cart(:,2));
legend("COM","Edge of Board")
xlabel("Position (x), meters")
ylabel("Position (y), meters")
title("Snowboarder on Ramp Position vs. Time")
hold off
%options2 = odeset('RelTol',reltol,'Events',event_stop_2)
figure(2)
plot(T2,body_angle_2)

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
%Stage 2 - Snowboarder in air's COM EOM - changing MOI
%z1 = x, z2 = y, z3 = xdot, z4 = ydot, z5 = t_b angle, z6 = t_bdot, z7=L

dz1dt = Z(3);
dz2dt = Z(4);
dz3dt = 0;
dz4dt = -g;
dz5dt = Z(6);
if T < (T1_end + MOI_change_time)
    dz6dt = 2*Ldot*Z(6)/Z(7);
else 
    dz6dt = 0;
end
dz7dt = Ldot;

dzdt = [dz1dt;dz2dt;dz3dt;dz4dt;dz5dt;dz6dt;dz7dt];

end
%{
Stage 3 EOM where MOI is always constant
function dzdt = eom3(T,Z)
    dzdt = [Z(3);Z(4);0;-g;Z(6);0];
end
%}    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventvalue,stopthecalc,eventdirection] = event_stop_1(T,Z)
    eventvalue      =  (Z(1)-((ramp_theta_cutoff + 90)*2*pi/360));%event at end of ramp defined by ramp angle cutoff
    stopthecalc     =  1;       %  Stop if event occurs
    eventdirection  = 1;       %  Detect only events with dydt>0
end
%{
function [eventval,stop,dir] = event_stop_2(T,Z)
    eventval = T - MOI_change_time;
    stop = 1;%stop if event occurs setting = true
    dir = 1;%event only occurs if dLdt > 0
end
%}
function [eventval,stop,dir] = event_stop_3(T,Z)
    eventval = Z(2); %event at y = 0
    stop = 1;%stop if event occurs setting = true
    dir = -1;%event only occurs if dydt < 0

end

end