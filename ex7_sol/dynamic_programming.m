clc; clear; close all;

Ts = 0.1;                   % time discretization grid constant
nSteps = 1;                 % number of integration steps during Ts
Q = diag([100,0.01]);       % state weighting matrix
R = 0.001;                  % control weighting matrix

% limits of state and control space
x1_max = 2*pi;
x1_min = -pi/2;
x2_max = 10;                % lower bound is -x2_max
u_max = 10;                 % lower bound is -u_max

% design LQR controller
% output written to lqr.mat
LQR_design(Q, R, Ts, nSteps);

N = 20;                     % control horizon
input.Ts = Ts;              
input.nSteps = nSteps; 
input.sens = 0;             % no sensitivities are needed here

N_x1 = 200;                 % number of discretization points for state x1
N_x2 = 40;                  % number of discretization points for state x2
N_u = 20;                   % number of discretization points for control u

% discretized state and control space
x1_values = linspace( x1_min, x1_max, N_x1);
x2_values = linspace(-x2_max, x2_max, N_x2);
u_values = linspace(-u_max, u_max, N_u);

% mesh for plotting
[X1, X2] = meshgrid(x1_values, x2_values);

% THESE VARIABLES ARE GENERATED AND STORED BY RUNNING LQR_design.m
load lqr.mat A B Q R K P

% initialize the cost-to-go function with the cost matrix associated with 
% LQR controller

LQR_cost = zeros(N_x1, N_x2);   % grid with cost of every state combination
LQR_u = zeros(N_x1, N_x2);      % LQR control for every state combination

% loop through all state combinations
for i = 1:N_x1
    for j = 1:N_x2
        x = [x1_values(i); x2_values(j)];
        % INSERT YOUR CODE HERE:
        LQR_cost(i,j) = x' * P * x;
        LQR_u(i,j) = -K*x;
    end
end

% terminal cost of OCP is same as LQR
J_cost = LQR_cost;

% compute cost-to-go of initial state via backward recursion
for k = N-1:-1:1
    
    % will contain optimal control input for every state
    u_map = NaN * ones(N_x1, N_x2);
    % initialization of cost-to-go function at current iteration
    % infinity until we know a better cost
    J_new = inf + J_cost;             

    % loop through all state combinations
    for i1 = 1:N_x1
        for i2 = 1:N_x2
%             cost_local = zeros(N_u,1);
            x_k = [x1_values(i1); x2_values(i2)];
            input.x = x_k;
            % loop through all control
            for j = 1:N_u
                u_k = u_values(j);
                input.u = u_k;
                % apply control
                output = RK4_integrator( @ode, input );
                
                % project on discretization grid
                % INSERT YOUR CODE HERE:
                i1_next = project(output.value(1), x1_values);
                i2_next = project(output.value(2), x2_values);
                
                % if not on grid
                if i1_next <= 0 || i1_next > N_x1 || i2_next <= 0 || i2_next > N_x2
                    cost = Inf;
                else
                    % INSERT YOUR CODE HERE:
                    % cost of this control at this state
                    cost = x_k' * Q * x_k + u_k' * R * u_k + J_cost(i1_next, i2_next) ;
                end
                
%                 cost_local(j) = cost;
                
                if cost < J_new(i1,i2)
%                   INSERT YOUR CODE HERE:
                    J_new(i1,i2) = cost;
                    u_map(i1,i2) = u_k;
                end
            end
        end
    end
    
    J_cost = J_new;
    
    figure(1);
    clf;
    subplot(211);
    surf(X1.', X2.', J_cost); hold all;
    plot3(X1.', X2.', LQR_cost,'--r');
    xlabel('x_1'); ylabel('x_2'); zlabel('J_{cost}');
    legend('DP', 'LQR')
    title(['Cost-to-go function for k = ' num2str(k)])
    drawnow
    
    subplot(212);
    surf(X1.', X2.', u_map); hold on;
    plot3(X1.', X2.', LQR_u,'--r');
    xlabel('x_1'); ylabel('x_2'); zlabel('u_{map}');
    zlim([-u_max; u_max])
    legend('DP', 'LQR')
    title(['Optimal feedback control for k = ' num2str(k)])
    view(150, 75);
    drawnow
    
end

save DP.mat x1_values x2_values u_values J_cost u_map Ts nSteps u_max

