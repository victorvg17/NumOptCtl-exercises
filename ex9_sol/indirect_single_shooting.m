close all;
clear;
clc;
%%
import casadi.*

T = 10;                             % integration time

% state variables
x = SX.sym('x', 2);                 % state
lam = SX.sym('lam', 2);             % costate
aug = [x; lam];                     % augmented state

% boundary values
x0bar = [0; 1];
lamTbar = [0;0];

% optimal control
uopt = max(-1, min(1, -lam(1) / 2));

% state dynamics
xdot = [(1 - x(2)^2) * x(1) - x(2) + uopt;
          x(1) ];

% costate dynamics
lamdot = [-2 * x(1) - lam(1) * (1 - x(2)^2) - lam(2);
          -2 * x(2) + 2 * lam(1) * x(1) * x(2) + lam(1)];

% augmented dynamics
augdot = [xdot; lamdot];

% integrator
% aim: create integrator function F(x_in, tf) that integrates the
% given dynamics starting at state x_in over a time interval tf and returns
% the resulting state.
tf = SX.sym('tf');                  % integration interval

% create struct with information needed to define the integrator
% the dynamics are scaled by tf.
% this is because CVODES by default integrates over the time interval [0,1].
% integration of the scaled dynamics from 0 to 1 is equivalent to
% integration of the original unscaled dynamics from 0 to tf
dae = struct('x', aug, 'p', tf, 'ode', tf * augdot);
% some options (precision, in the sense of absolut and relative error
% tolerance)
opts = struct('abstol', 1e-8, 'reltol', 1e-8);
% create the desired integrator function
F = integrator('F', 'cvodes', dae, opts);

%% build dummy NLP

% initial value of lambda (to be found)
lam0 = MX.sym('lam0', 2);

% compute augmented state at time T dependend on lam0
intout = F('x0', [x0bar; lam0], 'p', T);
augT = intout.xf;                      % augmented state at T 
lamT = intout.xf(3:4);                 % final value lambda(T)

% terminal condition
% this residual should be zero
g = lamT - lamTbar;
lbg = 0;
ubg = 0;

% NLP with dummy objective
nlp = struct('x', lam0, 'f', 0, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp);
sol = solver('lbg', lbg, 'ubg', ubg);

lam0opt = full(sol.x);

%% simulate solution

% integrate in N timesteps to get intermediate results
N = 100;
DT = T/N;

% integration loop
AUG = [x0bar; lam0opt];
for i = 1:N
    intres = F('x0', AUG(:, end), 'p', DT);
    AUG = [AUG, full(intres.xf)];
end

% split int state, costate, controls
xopt = AUG(1:2, :);
lamopt = AUG(3:4, :);
u_opt = max(-1, min(1, -lamopt(1,:) / 2));
% time grid
tvec = 0:DT:T;

figure(1); clf;
subplot(3,1,1); hold on;
plot(tvec, xopt(1,:) )
plot(tvec, xopt(2,:) )
legend('x_1', 'x_2')

subplot(3,1,2); hold on;
plot(tvec, lamopt(1,:) )
plot(tvec, lamopt(2,:) )
legend('\lambda_1', '\lambda_2')

subplot(3,1,3); hold on;
plot(tvec, u_opt)
legend('u')
xlabel('time t')
