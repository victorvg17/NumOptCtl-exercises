%% 2a) casadi expression for P
% Load CasADi
import casadi.*

v = MX.sym('v', 4);
pos_t = Function('pos_t', {v}, {ballistic_dynamics_RK4(v)});
v_start = [0, 0, 0, 0];
P = full(pos_t(v_start));

%% 2b) dynamic optimisation problem
% Create NLP:
% minimize  (X(1)-X(5))^2 + (X(2)-X(6))^2
% subject to X(2) >= 0
% subject to X(6) >= 0
% subject to X(3)^2 + X(3)^2 <= v_bar^2;
% subject to X(7)^2 + X(8)^2 <= v_bar^2;
v_bar = 15;
p0 = [0, 0, 10, 0];
p = SX.sym('p', 1, 4);
v = SX.sym('v', 1, 4);
k = [v, p];
f = full(pos_t(v));
g = [p(2);
     p(4);
     v(1)^2 + v(2)^2;
     v(3)^2 + v(4)^2];
 

% Create IPOPT solver object
nlp = struct('x', k, 'f', f, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp);

% Solve the NLP
res = solver('k0', [[0, 0, 0, 0],[0, 0, 10, 0]],...
             'lbx', -inf,...
             'ubx', inf,...
             'lbg', [0; 0; -inf; -inf],...
             'ubg', [inf; inf; v_bar^2; v_bar^2]);
         
% Print the solution
f_opt_ball = full(res.f);
k_opt_ball = full(res.k);
lam_x_opt_ball = full(res.lam_x);
lam_g_opt_ball = full(res.lam_g);