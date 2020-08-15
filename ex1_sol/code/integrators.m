clc; clear; close all;

dynamics = @(t,x) [x(2); -.2*x(2) - x(1)];

Ts = .5;
N = 19;
x0 = [1;0];
t = 0:Ts:Ts*N;

% 2a)
[t,X_ode45] = ode45( dynamics, t, x0);
X_ode45 = X_ode45';

figure(1); clf; hold on;
plot(t,X_ode45(1,:), 'kx-', 'DisplayName', 'x1 ode45')
plot(t,X_ode45(2,:), 'ko-', 'DisplayName', 'x2 ode45')
legend('Location','eastoutside');
xlabel('time t')
ylabel('state x')


%% 2b) euler

% remove formal time dependency of dynamics
dynamics = @(x) dynamics(0,x);

X_eu = x0;
h = .125;
t_eu = 0:h:Ts*N;
for k = 2:length(t_eu)
    X_eu(:,k) = X_eu(:,k-1) + h * dynamics(X_eu(:, k-1)) ;
end

plot(t_eu,X_eu(1,:), 'rx-', 'DisplayName', 'x1 euler')
plot(t_eu,X_eu(2,:), 'ro-', 'DisplayName', 'x2 euler')

%% 2b) rk4
X_rk4 = x0;
for k = 2:length(t)
    X_rk4(:,k) = rk4step(dynamics, Ts,  X_rk4(:,k-1));
end   

plot(t,X_rk4(1,:), 'bx-', 'DisplayName', 'x1 rk4')
plot(t,X_rk4(2,:), 'bo-', 'DisplayName', 'x2 rk4')

%% 2c) casadi 
import casadi.*
x = MX.sym('x',2);
dt = MX.sym('x',1);
x_next = Function('x_next',{x,dt}, {rk4step(dynamics, dt, x)});

X_ca = x0;
for k = 2:length(t)
    X_ca(:,k) = full(x_next( X_ca(:,k-1), Ts)) ;
end 

plot(t,X_rk4(1,:), 'gx-', 'DisplayName', 'x1 casadi')
plot(t,X_rk4(2,:), 'go-', 'DisplayName', 'x2 casadi')
