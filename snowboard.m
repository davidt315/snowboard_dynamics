%
function snowboard
clear all
close all


g = 9.81;                   % gravitational acceleration in m/s^2
m = 50;                     % mass in kg
r = 20;                     % ramp radius in m
ramp_theta_cutoff = 40;     % where ramp ends, degrees


%Model the snowboarder as a cylinder rotating around an axis through the
%"long" side, not the circular faces. The snowboarder has an MOI of 1/12
%ml^2, where l = the distance from their COM to the board. We change the
%length of the cylinder at a constant rate.
COM_board_dist = 0.5*1.754;    % distance btwn snowboard/ground and person's COM: Avg height*0.5
COM_board_end = 0.8382;         % distance btwn snowboard/ground and person's COM, crouched
MOI_change_time = 1;           % time taken for moment of inertia to change
Ldot = (COM_board_end - 2*COM_board_dist)/MOI_change_time;
r_COM = r - COM_board_dist;

% initial values phase 1
z0_1=[0, 0];

% Time and ode
t_span = [0, 12]; 

reltol = 1.0e-9;
options1= odeset('RelTol', reltol,'Events',@event_stop_1);

[T1,Z1] = ode45(@eom1, t_span, z0_1, options1);

T1_end = T1(end);

vphase1com = r_COM.*Z1(:,2);
vphase1board = r.*Z1(:,2);

%Stage 2 setup - defining IC's etc
x0 = r + r_COM*sind(ramp_theta_cutoff);     %considering left edge of ramp where x = 0
y0 = r - r_COM*cosd(ramp_theta_cutoff);     %considering lowest point of ramp where y=0
v0_mag = r_COM*(Z1(length(Z1),2));
vx0 = v0_mag*cos(Z1(length(Z1),1) - pi/2);
vy0 = v0_mag*sin(Z1(length(Z1),1) - pi/2);

%angle of snowboarder's body - COM relative to edge
t_b0 = Z1(length(Z1),1);
t_bdot0 = Z1(length(Z1),2);
z0_2 = [x0,y0,vx0,vy0,t_b0,t_bdot0,COM_board_dist];
t_span2 = [T1_end, 10];
options2 = odeset('RelTol',reltol,'Events',@event_stop_3);

[T2,Z2] = ode45(@eom2,t_span2,z0_2,options2);

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

% Full trajectory of snowboarder
hold on
plot(COMcart(:,1),COMcart(:,2));
plot(Edge_cart(:,1),Edge_cart(:,2));
plot([COMramp_cart(length(COMramp_cart),1) Edge_ramp_cart(length(Edge_ramp_cart),1)], ...
    [COMramp_cart(length(COMramp_cart),2) Edge_ramp_cart(length(Edge_ramp_cart),2)], '-black');
legend("COM","Edge of Board","End of Ramp")
xlabel("Position (x), meters")
ylabel("Position (y), meters")
title("Snowboarder on Ramp Position vs. Time")
hold off


% Angular Velocity and Accelerations
% xbody = Z2(:,1);
% ybody = Z2(:,2);
% xramp = Z2(:,1) - Z2(:,7);
% yramp = Z2(:,2) - Z2(:,7);
vxbody = Z2(:,3);
vybody = Z2(:,4);
vbody = sqrt(vxbody.^2 + vybody.^2);
theta = Z2(:,5);
thetadot = Z2(:,6);
Ltime = Z2(:,7);
vxboard = vxbody + Ltime.*thetadot.*sin(theta);
vyboard = vybody - Ltime.*thetadot.*cos(theta);
vboard = sqrt(vxboard.^2 + vyboard.^2);

merged_vbody = [vphase1com; vbody];
merged_vboard = [vphase1board; vboard];
merged_thetadot = [Z1(:,2); thetadot];


figure()
subplot(2,1,1)
plot(full_time, merged_vbody)
hold
plot(full_time, merged_vboard)
hold
xline(T1_end)
hold
ylabel('Linear Velocity')
title('Linear and Angular Velocities')
legend('COM Velocity', 'Board Velocity', 'Ramp Cutoff', 'Location', 'best')
subplot(2,1,2)
plot(full_time, merged_thetadot)
hold
xline(T1_end)
hold
ylabel('Angular Velocity')
xlabel('Time (s)')
legend('COM Angular Velocity', 'Ramp Cutoff', 'Location', 'best')
figure(2)


% Energy
PEbody = m.*g.*COMcart(:,2);
angvel_e = [0.*T1; thetadot.^2.*((1/12).*m.*Z2(:,7))];
KEbody = m.*merged_vbody.^2./2 + angvel_e;
TEbody = PEbody + KEbody;

figure();
hold on;
plot(full_time, PEbody);
plot(full_time, KEbody);
plot(full_time, TEbody);
title('COM Energy Over Time')
xlabel('Time (s)')
ylabel('Energy (J)')
legend("Potential Energy", "Kinetic Energy", "Total Energy", 'Location', 'best');
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = eom1 (T,Z)
%Stage 1 (Snowboarder on ramp EOM)
% DAE form
% eom + constraint eqn (written wrt state variables) 
%z1 = theta, z2 = thetadot
F_d = 0;%adding a drag force = 0 in case we want to implement this later

dz1dt = Z(2);
dz2dt = g*cos(Z(1))/r_COM - F_d/m;
dzdt = [dz1dt;dz2dt];

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
    dz6dt = -2*Ldot*Z(6)/Z(7);
    dz7dt = Ldot;
else   
    dz6dt = 0;
    dz7dt = 0;
end


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