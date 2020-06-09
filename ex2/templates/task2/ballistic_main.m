clc; clear; close all;

import casadi.*
%parameters
vbar = 15;
w0 = zeros(4, 1);

%decision variables
v1 = SX.sym('v1', 2);
v2 = SX.sym('v2', 2);
w = [v1; v2];

%objective
P = ballistic_dynamics_RK4(w)';
p1 = P(1:2);
p2 = P(3:4);
dp = p1 - p2;
f = dp'*dp;

%% constraints
g = [];
lbg = [];
ubg = [];

% p1z >= 0
g = [g; p1(2)];
lbg = [lbg; 0];
ubg = [ubg; inf];

% p2z >= 0
g = [g; p2(2)];
lbg = [lbg; 0];
ubg = [ubg; inf];

% v1^2 <= vbar^2
g = [g; v1' * v1];
lbg = [lbg; -inf];
ubg = [ubg; vbar^2];

% v2^2 <= vbar^2
g = [g; v2' * v2];
lbg = [lbg; -inf];
ubg = [ubg; vbar^2];

nlp = struct('x', w, 'f', f', 'g', g);

%% Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);

% Solve the NLP
res = solver('x0' , w0,...            % solution guess
             'lbx', -inf,...          % lower bound on x
             'ubx', inf,...           % upper bound on x
             'lbg', lbg,...           % lower bound on g
             'ubg', ubg);             % upper bound on g

% simulate optimal solution
wsol = full(res.x);
T = 0.5;
M = 100;
DT = T/M;
X0 = [0, 0, wsol(1), wsol(2), 10, 0, wsol(3), wsol(4)]';
X = X0;

% RK4 integrator
for j=1:M
    X(:,j+1) = rk4step(@(t,x) ode(x), DT, X(:,j), 0 ) ;
end

figure(1); clf; hold on;
plot(X(1,:), X(2,:))
plot(X(5,:), X(6,:))