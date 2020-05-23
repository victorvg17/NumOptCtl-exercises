function [t_next, x_next] = rk_integrator(t_span, x0, h)
x_next = [];
t_next = [];
t_l = t_span(1);
t_h = t_span(2);

x_n = [x0(1), x0(2)];
t_n = t_l;
while t_n <= t_h
    k1 = myderi_fun(x_n);
    k2 = myderi_fun(x_n + (h/2)*k1);
    k3 = myderi_fun(x_n + (h/2)*k2);
    k4 = myderi_fun(x_n + h*k3);
    
    x_prime = x_n +(h/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    x_next = [x_next; [x_prime(1), x_prime(2)]];
    t_next = [t_next; t_n];
    
    x_n = x_prime;
    t_n = t_n + h;
end