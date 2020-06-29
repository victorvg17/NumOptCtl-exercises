function dx = ode(x,u)

    dx = [  x(2); ...
            sin(x(1)) + u(1)];

end

