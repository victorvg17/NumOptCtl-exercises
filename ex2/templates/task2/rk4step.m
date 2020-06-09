function x_next = rk4step(ode, h, x0, t0 )

x0 = x0(:) ;

k1 = ode(t0, x0);
k2 = ode(t0 + h/2, x0 + h/2 * k1);
k3 = ode(t0 + h/2, x0 + h/2 * k2);
k4 = ode(t0 + h, x0 + h * k3);

x_next = x0 + h/ 6 * (k1 + 2*k2 + 2*k3 + k4);
