%
function snowboard
clear all
close all


g = 9.81;                   % gravitational acceleration in m/s^2
m = 50;                     % mass in kg
r = 20;                     % ramp radius in m
ramp_theta_cutoff = 15;     % degrees where ramp ends

% initial values phase 1
z0=[0, 0, 0, 0];

% Time and ode
t_span = [0, 8]; 

reltol = 1.0e-6;
options= odeset('RelTol', reltol);

[T,Z] = ode45(@eom1, t_span, z0, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = eom1 (T,Z)
% DAE form
% eom + constraint eqn (written wrt state variables) 
% z1 = theta1, z2 = theta2, z3 = theta1dot, z4 = theta2dot, z5 = T1, z6 =T2

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

%
dzdt = Mat\vec;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end