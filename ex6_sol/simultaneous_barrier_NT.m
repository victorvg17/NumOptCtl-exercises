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

% collect all decision variables in w
g = [];

w = {};
w0 = 0.1;
L = 0;
roh = .01;

X_k = x0bar;

for i = 1:N
    U_k = MX.sym(['U_',num2str(i)],   nu, 1);
    L = L + sum(X_k.^2);
    L = L + 2 * sum(U_k.^2);
    L = L - roh * (log(-U_k + 1) + log(U_k + 1)); 

    X_next = F(X_k, U_k);
    
    X_k = MX.sym(['X_',num2str(i+1)], nx, 1);
    w = {w{:}, U_k};
    w = {w{:}, X_k};
    g = [g; X_next - X_k];
end

L = L + 10 * sum(X_k.^2);

z = {};
Lagr_expr = L;
for i = 1:N
    lam_k = MX.sym(['lam_',num2str(i)], nx, 1);     % lagrangian multipliers
    Lagr_expr = Lagr_expr + lam_k' * g(2*i-1:2*i);  % Lagrangian contrib
    z = {z{:}, w{2*i-1}};
    z = {z{:}, lam_k};
    z = {z{:}, w{2*i}};
end
z = vertcat(z{:});

Lagr = Function('Lagr', {z}, {Lagr_expr});
gradL = Function('gradL', {z}, {jacobian(Lagr_expr, z)'});
hessL = Function('hessL', {z}, {hessian(Lagr_expr, z)});

% newton type solver
max_it = 100;
rhs_tol = 1e-6;
iter = .1 * ones(size(z));

for i=2:max_it
    
    rhs = full(gradL(iter(:,i-1)));
    KKT = full(hessL(iter(:,i-1)));
    if norm(rhs) < rhs_tol
        fprintf('converged after %d iterations\n', i)
        break
    end
    wk = iter(:,i-1);
    step =  - KKT \ rhs;
    % line search
    alpha = 1;
    beta = .8;
    while true
        candidate = wk + alpha * step;
        U = candidate(1:2*nx+nu:end);
        if all(U < 1) && all(U > -1)
            break
        end
        alpha = beta * alpha;
    end
    iter(:,i) = candidate;
end
if i == max_it
    display('maximum iterations reached')
end

%% visualize solution

w_opt = iter(:, end);
U_opt = w_opt(1:2*nx+nu:end);
X_opt = [w_opt(4:2*nx+nu:end)';...
         w_opt(5:2*nx+nu:end)'];

X_opt = [x0bar, X_opt];

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