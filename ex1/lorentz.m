clc; clear; close all;
rho = 28;
beta = 8/3;
sigma = 10;

x0 = [1; 0; 0];
h = 0.01;
t = 0:h:100;

ode = @(x) [sigma*(x(2) - x(1));
            x(1)*(rho - x(3)) - x(2);
            x(1)*x(2) - beta*x(3)];
X = x0;
for k = 2:length(t)
    X(:, k) = rk4step(ode, h, X(:, k-1));
end


figure(1); clf; hold on;
plot3(X(1, :), X(2, :), X(3, :));
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('lorenz attractor');
