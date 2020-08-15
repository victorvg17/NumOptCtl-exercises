clc; clear; close all;

rho = 28;
sigma = 10;
beta = 8/3;

ode = @(x) [  sigma * (x(2) - x(1)); 
                x(1) * (rho - x(3)) - x(2);
                x(1) * x(2) - beta * x(3) ];
            
x0 = [1; 0; 0];
h = .01;
t = 0:h:100;

X = x0;
for k = 2:length(t)
    X(:,k) = rk4step(ode, h,  X(:,k-1));
end   

figure(1); clf; hold on;

plot3( X(1,:), X(2,:), X(3,:) )
