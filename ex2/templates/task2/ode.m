function [ dx ] = ode( x )

% parameters
d1 = 0.1;
d2 = 0.5;


g = 9.81;
w = 1;


% x_dot = f(x)
dx = [  x(3); ...
        x(4); ...
        (-(x(3)-w)*sqrt((x(3)-w)^2 + x(4)^2)*d1); ...
        (-x(4)*sqrt((x(3)-w)^2 + x(4)^2)*d1 - g); ...
        x(7); ...
        x(8); ...
        (-(x(7)-w)*sqrt((x(7)-w)^2 + x(8)^2)*d2); ...
        (-x(8)*sqrt((x(7)-w)^2 + x(8)^2)*d2- g); ];


end

