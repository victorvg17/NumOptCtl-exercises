function dxdt = derivative_fun(t, x)
dxdt = [x(2); -0.2*x(2)-x(1)];