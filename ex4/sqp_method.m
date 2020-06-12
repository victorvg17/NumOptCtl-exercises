%% Implementation of a Gauss-Neston SQP solver usign CasADi

close all
clear variables
clc

% import CasADi
import casadi.*

nv = 2;
x = MX.sym('x',nv);

% Define objective function
% Insert your code here:
% F = ...;
% Jf = ...;

% Define residuals (for Gauss-Newton Hessian approximation)
% Insert your code here:
% R = ...;
% Jr = ...;

% Define equalities 
% Insert your code here:
% G = ...;
% Jg = ...;

% Define inequalities 
% Insert your code here:
% H = ...;
% Jh = ...;

% linearization point (parameter)
xk = MX.sym('xk',nv);

% decision variable
delta_x = MX.sym('delta_x',nv);

% Define linearized constraints 
% Insert your code here:
% linearized equality constraints
% g_l = ...;
% linearized inequality constraints
% h_l = ...;

% Gauss-Newton Hessian approximation
% Bk = ...
% define the quadratic Gauss-Newton objective function
% f_gn = ...;

% Allocate QP solver
qp = struct('x',delta_x, 'f',f_gn,'g',[g_l;h_l],'p',xk);
solver = qpsol('solver', 'qpoases', qp);

% SQP solver
max_it = 100;
xk = [-2 4]'; % Initial guess
iter = zeros(nv,max_it);
iter(:,1) = xk;

lbg =  [0;-Inf];
ubg =  [0;0];

    
for i=2:max_it    
    
    % Solve the QP
    sol = solver('lbg', lbg, 'ubg', ubg , 'p', xk);
    step = full(sol.x);
    if norm(step) < 1.0e-16
        break;
    end
    
    % Line-search 
    t = 1;
    kappa = 0.7;
    alpha = 1.1;
    out = F(iter(:,i-1));
    prev_cost = full(out);
    next_cost = Inf;
    while(next_cost > alpha*prev_cost)
        trial = iter(:,i-1) + t*step;
        out = F(trial);
        next_cost = full(out);
        t = t*kappa;
    end
    iter(:,i) = iter(:,i-1) + t*step;
    xk = iter(:,i);
    
end
iter = iter(:,1:i-1);

figure()
% Plot feasible set, and iterations
pl = ezplot('sin(x) - y^2');
set(pl,'Color','red');
hold all
pl= ezplot('x^2 + y^2 - 4');
set(pl,'Color','blue');
ezcontour('(x- 4)^2 + (y - 4)^2')
plot(iter(1,:),iter(2,:),'--','Color','black')
plot(iter(1,:),iter(2,:),'o','Color','black')
title('Iterations in the primal space')
grid on

figure()
plot(iter(1:2,:)')
grid on
xlabel('iterations')
ylabel('primal solution')
grid on