% code snippet for how to build a multiple shooting OCP in casadi
% adapted from direct_multiple_shooting.m from the casadi example pack
% for a complete example, but with lots of details you will not need, see
% https://github.com/casadi/casadi/releases/download/3.4.5/casadi-example_pack-v3.4.5.zip


%% declare all functions and parameters you need to formulate the NLP
clc; clear; close all;
N = 50;        % control horizon
nx = 2;        % state dimension
nu = 1;        % control dimension
h = 0.1;
import casadi.*
% xk = MX.sym('xk', nx, 1);
% uk = MX.sym('uk', nu, 1);
dynamics = @ (x, u)[x(2); sin(x(1) + u)];  % ode of the systems
F = @(x, u) rk4step(x, u, dynamics, h); %integrator from x_k, u_k to x_k+1
...

%% formulate and solve the NLP
% Start with an empty NLP
w = {};             % decision variables
J = 0;              % cost
g = {};             % constraints

% elimination of initial state -> x0 is not a decision variable
x0bar = [pi; 0];
xk = x0bar;

% build decision variables, objective, and constraint
for k = 0:N-1
    % New NLP variable for the control u_k
    uk = MX.sym(['u_', num2str(k)], nu);
    
    % collect in w
    w = {w{:}, uk};            

    % Integrate till the end of the current interval
    % something something xk uk
%     xnext = rk4step(xk, uk, dynamics, h); 
    xnext = F(xk, uk); 
    
    % contribution of stage cost to total cost
    J = J + (xk'*xk + 2*uk^2);

    % New NLP variable for state at end of interval
    xk = MX.sym(['x_', num2str(k+1)], nx);
    
    % collect in w
    w = {w{:}, xk};
    
    % Add dynamic constraint function,
    % constraining xk to integration result
    g = {g{:},  xk-xnext};

end

% contribution of terminal cost
J = J + 10*(xnext'*xnext);

%% structure of hessian of lagrangian
% NOT A NECESSARY PART TO BUILD NLP
% % Lagrangian
% L = J;          % start with contribution of objective function
% 
% % collect all variables of lagrangian in this (x_k, u_k, lambda_k) in the
% % order defined on the sheet
% z = {};
% for k = 1:N
%     % New variable for multiplier lambda_k of k-th dynamic constraint
%     lam_k = MX.sym(['lam_', num2str(k)], nx);
%     
%     % contribution of constraint k to lagrangian
%     L = L + ...something something g{k}, lam_k
%         
%     % collect variables in correct order
%     z = {z{:}, w{?}};   % u_k-1
%     z = {z{:}, lam_k};
%     z = {z{:}, w{?}};   % x_k
% end
% z = vertcat(z{:});          % transform to column vector
% 
% % Hessian of Lagrangian
% Hess_L = Function('Hess_L', {z}, ...)
% ... evaluate at some value and spy
% END NOT NECESSARY PART

% Create an NLP solver
% vertcat({w}) will put all elements of w into a column vector
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% either build lbg and ubg along with g, but as vector, like
% lbg = [];
% lbg = [lbg; ...];
% lbg = [lbg; ...];
% ...
% or just use g in standard form g(x) = 0 and use
lbg = 0; ubg = 0;

% similar for w0
% for complicated initial guess build along side w, otherwise just
% put it to some value here.
w0 = zeros(N*(nu+nx), 1);

% Solve the NLP
sol = solver('x0', w0, 'lbg', lbg, 'ubg', ubg);

% obtain solution
w_opt = full(sol.x);
u_opt = w_opt(1:nx+nu:end);
x1_opt = w_opt(2:nx+nu:end);
x2_opt = w_opt(3:nx+nu:end);

%%  visualize solution
t_oc = 0:h:h*(N-1);
figure(1); clf;
plot(t_oc, u_opt, '-+', 'DisplayName', 'control trajectory');
legend('Location', 'east outside');
xlabel('time t');
ylabel('control u');

figure(2);
plot(x1_opt, x2_opt, '-+', 'DisplayName', 'state trajectory');
legend('Location', 'east outside');
xlabel('x1');
ylabel('x2');