function dxdt = myderi_fun(x)
dxdt = [x(2), -0.2*x(2)-x(1)];