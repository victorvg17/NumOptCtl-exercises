clc; clear; close all;
dynamics = @ (t, x) [x(2); -0.2*x(2) - x(1)];

x0 = [1; 0];
Ts = 0.5;
N = 19;
t = 0:Ts:Ts*N;

%% 2a)
[t, X_ode45] = ode45(dynamics, t, x0);
X_ode45  = X_ode45';

figure(1); clf; hold on;
plot(t, X_ode45(1, :), '-+', 'DisplayName', 'x1 ode45');
plot(t, X_ode45(2, :), '-o', 'DisplayName', 'x2 ode45');
legend('Location', 'east outside');
xlabel('time t');
ylabel('state x');

%% 2b) explicit euler
% remove time dependancy
dynamics = @(x) dynamics(0, x);

X_eu = x0;
h = 0.125;
t_eu = 0:h:Ts*N;

for k = 2:length(t_eu)
    X_eu(:, k) = X_eu(:, k-1) + h*dynamics(X_eu(:, k-1));
end
plot(t_eu, X_eu(1, :), '-+', 'DisplayName', 'x1 euler');
plot(t_eu, X_eu(2, :), '-o', 'DisplayName', 'x2 euler');

%% 2b) rk4
X_rk4 = x0;
h = 0.125;
t_rk4 = 0:h:Ts*N;

for k = 2:length(t_rk4)
    X_rk4(:, k) = rk4step(dynamics, h, X_rk4(:, k-1));
end

plot(t_rk4, X_rk4(1, :), '-+', 'DisplayName', 'x1 rk4');
plot(t_rk4, X_rk4(2, :), '-o', 'DisplayName', 'x2 rk4');

%% 2c) casadi
import casadi.*
x = MX.sym('x', 2);
dt = MX.sym('t');
x_next = Function('x_next', {x, dt}, {rk4step(dynamics, dt, x)});

X_ca = x0;
for k = 2:length(t)
    X_ca(:, k) = full(x_next(X_ca(:, k-1), Ts));
end

plot(t, X_ca(1, :), '-+', 'DisplayName', 'x1 ca');
plot(t, X_ca(2, :), '-o', 'DisplayName', 'x2 ca');
