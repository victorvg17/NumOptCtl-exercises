clear;  clc; close all;
import casadi.*

% set to 1 / true to solve NLP for many alpha and plot the gradients
% of the constraints and the lagrange multipliers
PLOT_GRADIENTS = 1;

if PLOT_GRADIENTS
    n_points = 100;         % num of alpha values for which NLP is solved
    alpha_min = -1;
    alpha_max = 1;
    alpha_v = linspace(alpha_min,alpha_max,n_points);
else
    n_points = 1;
    alpha_v = 0;            % solve NLP for this alpha only
end

T = 0.5;        % time horizon
N = 100;        % discretization steps
p0 = [0;0];     % Initial position of the ball
v_max = 15;     % max velocity

alpha = MX.sym('alpha',1,1);    % variable for alpha parameter

% TODO: complete file ballistic_dynamics.m
x = MX.sym('x',4,1);            % state variable
xdot = ballistic_dynamics(x);   % Model equations
f = Function('f', {x}, {xdot}); % ode function 

%% Formulate discrete time dynamics
DT = T/N;

% One step Runge-Kutta4 integrator
k1 = f(x);
k2 = f(x + DT/2 * k1);
k3 = f(x + DT/2 * k2);
k4 = f(x + DT * k3);
xnext = x + DT/6 * (k1 + 2*k2 + 2*k3 + k4);
rk4step = Function('rk4step', {x}, {xnext});

%% Build final state as expression of initial velocity

% initial position is fixed, initial velocity the decision variable
v0 = MX.sym('v0', 2, 1);        % initial velocity 
X = [p0; v0];                   % initial state

% repeated integration
for j=1:N
    X = rk4step(X);
end

% final state as function of initial velocity
F = Function('F', {v0}, {X}); 

%% NLP Formulation

% decision variables
w = v0;
w0 = [0.001; 0.001];
lbw = -inf;
ubw = inf;

% final state as expression of v0
xf = F(v0);

% cost
J = -xf(1);

% constraints
g = [];
lbg = [];
ubg = [];

% first constraint (equation (4b))
g = [ g; -xf(2) ];
lbg = [lbg; -inf];
ubg = [ubg; 0];

% second constraint (equation (4c))
g = [g; -alpha(xf(1) - 10) - xf(2)];
lbg = [lbg; -inf];
ubg = [ubg; 0];

% third constraint  (equation (4d))
g = [g; v0(1)^2 + v0(2)^2];
lbg = [lbg; -inf];
ubg = [ubg; v_max^2];

% Insert you code here (Jacobian of the constraints)
Jg = Function('Jg', {w, alpha}, {jacobian(g, w)});

% Create an NLP solver
% pass alpha as a parameter ('p')
prob = struct('f', J, 'x', w, 'g', g ,'p',alpha);
solver = nlpsol('solver', 'ipopt', prob);
% plotting preparations
% constrained norm of v0 constr
theta = linspace(0, 2*pi, 100);
Vcirc = v_max * [sin(theta); cos(theta)];

% ground constraint
vv = linspace(-15, 15, 20);                % plot constraints in this range
[V1, V2] = meshgrid(vv,vv);

G_ground = zeros(size(V1));
for i = 1:size(V1,1)
    for j = 1:size(V1,2)
        xxx =  F( [V1(i,j) ; V2(i, j)] );
        G_ground(i,j) = full(-xxx(2));
    end
end

lambdas = zeros(3,n_points);
for i = 1:n_points
    
    alpha_val = alpha_v(i);
    
    % Insert you code here (solve the NLP)
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
        'lbg', lbg, 'ubg', ubg, 'p', alpha_val);
    
    w_opt = full(sol.x);
    Jg_eval = full(Jg(w_opt,alpha_val));
    
    if PLOT_GRADIENTS
         % legend will depend on the order in which you defined the constraints
        figure(1)
        % plot the normalized gradients
        subplot(211);
        plot([0, Jg_eval(1,1).'/norm(Jg_eval(1,:))], [0, Jg_eval(1,2).'/norm(Jg_eval(1,:))]);
        hold on     % so the plot commands don't overwrite each other
        plot([0, Jg_eval(2,1).'/norm(Jg_eval(2,:))], [0, Jg_eval(2,2).'/norm(Jg_eval(2,:))]);
        plot([0, Jg_eval(3,1).'/norm(Jg_eval(3,:))], [0, Jg_eval(3,2).'/norm(Jg_eval(3,:))]);
        hold off    % so next plot() ins this subfig  will overwrite the old ones
        title('constraint Jacobians (normalized)')
        xlabel('y')
        ylabel('z')
        grid on
        xlim([-1,1])
        ylim([-1,1])
        legend('constr1','constr2','constr3','Location','eastoutside')
        
        % plot the lagrange multipliers
        subplot(212)
        lambdas(:,i) = full(sol.lam_g);
        semilogy(alpha_v(1:i) ,abs(lambdas(:,1:i).'));
        title('Lagrange multipliers')
        xlabel('\alpha')
        ylabel('\lambda')
        grid on
        legend('constr1','constr2','constr3','Location','eastoutside')
    end
end
   
if ~PLOT_GRADIENTS
    % compute trajectory
    X_traj = zeros(4, N+1);
    X_traj(:,1) = [p0; w_opt(1); w_opt(2);];    % initial state
    % integrate
    for i = 2:N+1
        X_traj(:,i) = full(rk4step(X_traj(:,i-1)));
    end
    
    plot(X_traj(1,:), X_traj(2,:))
    hold all
    xlabel('y')
    ylabel('z')
    grid on
end