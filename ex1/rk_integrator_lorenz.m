function [t_next, x_next] = rk_integrator_lorenz(t_span, x0, h, param)
x_next = [];
t_next = [];
t_l = t_span(1);
t_h = t_span(2);

x_n = [x0(1), x0(2), x0(3)];
t_n = t_l;
while t_n <= t_h
    k1 = lorenz_attractor_deri(x_n, param);
    k2 = lorenz_attractor_deri(x_n + (h/2)*k1, param);
    k3 = lorenz_attractor_deri(x_n + (h/2)*k2, param);
    k4 = lorenz_attractor_deri(x_n + h*k3, param);
    
    x_prime = x_n +(h/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    x_next = [x_next; [x_prime(1), x_prime(2), x_prime(3)]];
    t_next = [t_next; t_n];
    
    x_n = x_prime;
    t_n = t_n + h;
end