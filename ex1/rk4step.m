function x_next = rk4step(ode, h, x)
k1 = ode(x);
k2 = ode(x + k1*h/2);
k3 = ode(x + k2*h/2);
k4 = ode(x + k3*h);
x_next = x + h/6*(k1 + 2*k2 + 2*k3 + k4);