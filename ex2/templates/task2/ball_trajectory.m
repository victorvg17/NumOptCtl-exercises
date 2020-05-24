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

x = SX.sym('x', 1, 8);
x_init = [0, 0, 0, 0, 10, 0, 0, 0];
f = (x(1) - x(5))^2 + (x(2) - x(6))^2;
curr_pos = [x(1), x(2), x(5), x(6)];
curr_vel = [x(3), x(4), x(7), x(8)];
g = full(pos_t(curr_vel)) - curr_pos;

h = [-x(2); 
     -x(6);
      x(3)^2 + x(3)^2 - v_bar^2;
      x(7)^2 + x(8)^2 - v_bar^2];

nlp = struct('x', x, 'f', f, 'g', g, 'h', h);

% Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);

% Solve the NLP
res = solver('x0' , [0, 0, 0, 0, 10, 0, 0, 0],... % solution guess
             'lbx', -inf,...           % lower bound on x
             'ubx',  inf,...           % upper bound on x
             'lbg',    0,...           % lower bound on g
             'ubg',    0);             % upper bound on g

         % Print the solution
f_opt = full(res.f)          
x_opt = full(res.x)          
lam_x_opt = full(res.lam_x) 
lam_g_opt = full(res.lam_g)