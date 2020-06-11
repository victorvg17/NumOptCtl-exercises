function [ xdot ] = ballistic_dynamics( x ) 

    % ode parameters
    w = 2;          % wind speed
    d = .1;         % drag coefficient
    g = 9.81;       % gravity
    v = x(3:4);
    dnorm = norm(v - [w;0]);

    % Insert your code defining the ode here
    xdot = [v(1); v(2); ...
            (v(1) - w)*dnorm*d;
            v(2)*dnorm*d - g];

end

