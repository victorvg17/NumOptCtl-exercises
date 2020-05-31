%% 2a) differential equation for the system
%%x = [Py, Pz, Vy, Vz]' 
function [ dx ] = ode( x )

% parameters
d1 = 0.1;
g = 9.81;
w = 1;
M = 1;

% x_dot = [Py_dot, Pz_dot, Vy_dot, Vz_dot]'
dx = [  x(1); ...
        x(2); ...
        (-(x(1)-w)*sqrt((x(1)-w)^2 + x(2)^2)*d1); ...
        (-x(2)*sqrt((x(1)-w)^2 + x(2)^2)*d1 - g);];

end

