function dxdt = lorenz_attractor_deri(x, param)
dxdt = [param.sigma*(x(2) - x(1)), x(1)*(param.rho - x(3)) - x(2), x(1)*x(2) - param.beta*x(3)];