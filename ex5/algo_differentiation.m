close all
clear variables
clc

% import CasADi
import casadi.*

nv = 1;
x0 = 0.5;
N = 50;
h = 0.1;
uk = MX.sym('uk', N);

% phi = zeros(N+1);
phi = MX.sym('phi', N);
phi(1) = ((1-x0)*x0 + uk(1) )*h + x0;
for i = 2:N
   phi(i) = ((1-phi(i-1))*phi(i-1) + uk(i) )*h + phi(i-1);
end    