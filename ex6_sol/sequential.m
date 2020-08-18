close all; clear; clc;
import casadi.*;

%% initial

% parameters
nx = 2;         % state dimension
nu = 1;         % control dimension
N = 50;         % horizon length
DT = .1;        % discretization time step
N_rk4 = 10;     % number of rk4 steps per discretiazion time step
x0bar = [pi; 0];   % initial state

% build integrator
% dynamics
dynamics = @(x,u) [x(2); sin(x(1)) + u];
h = DT / N_rk4;             % integration step
x = MX.sym('x',nx,1);
u = MX.sym('u',nu,1);

x_next = x;
for i = 1:N_rk4
    x_next = rk4step(x_next, u, dynamics, h);
end

% integrator
F = Function('F',{x, u}, {x_next});

%% NLP formulation
U = MX.sym('U', N);         % vector of all controls
U0 = .1 * ones(N, 1);         % intial guess U
% build vector of all states dependend on U
X = F(x0bar, U(1));
for i = 1:N-1
    X = [X, F(X(:, i), U(i+1))];
end

% cost
L = sum(sum(X(:,1:end-1).^2)) + 2 * sum(U.^2);          % stage costs
L = L + 10 * sum(X(:,end).^2);                          % terminal cost


hessL = Function('hessL', {U}, {hessian(L, U)});

figure(1)
spy(full(hessL(U0)));

% create nlp solver
nlp = struct('x', U, 'f', L);
solver = nlpsol('solver','ipopt', nlp);

% solve nlp
sol = solver('x0', U0);

%% visualize solution

U_opt = sol.x;

FX = Function('FX', {U}, {X});
X_opt = [x0bar, full(FX(U_opt))];

figure(2); clf;
subplot(2,1,1); hold on;
plot(0:N, X_opt(1,:))
plot(0:N, X_opt(2,:))
title('state trajectory')
legend('\phi', '\omega')

subplot(2,1,2); hold on;
stairs(0:N, [full(U_opt); nan])
title('control trajectory')
legend('\tau')
xlabel('discrete time k')
%% animate

animatePendulum(X_opt(1,:), 0.05, 'pendulum.gif')
