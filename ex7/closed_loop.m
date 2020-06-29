clc; clear; close all;

input.sens = 0; % no sensitivities are needed here
x0 = [pi; 0];   % intial state

% THESE VARIABLES ARE GENERATED AND STORED BY RUNNING LQR_design.m and 
% dynamic programming.m
load lqr.mat A B Q R K P
load DP.mat x1_values x2_values u_values J_cost u_map Ts nSteps u_max

input.Ts = Ts;
 % do twice as many integration steps as the controller
input.nSteps = nSteps * 2; 

% Closed-loop simulation
time = 0;
Tf = 5;
state_LQR = x0;
state_DP = x0;
us_DP = [];
us_LQR = [];
cost_DP = 0;
cost_LQR = 0;
iter = 0;
figure(1);
u_DP_old = 0;

while time(end) < Tf
    
    % optimal feedback law from LQR
    % INSERT YOUR CODE HERE:
    % don't forget to clip to the contraints
    u_LQR = ...;
    us_LQR = [us_LQR u_LQR]; % remember control trajectory
    % add incurred cost
    cost_LQR = [cost_LQR cost_LQR(end)+u_LQR.'*R*u_LQR+state_LQR(:,end).'*Q*state_LQR(:,end)];
    
    % apply LQR control to nonlinear system
    input.x = state_LQR(:,end);
    input.u = u_LQR;
    output = RK4_integrator( @ode, input );
    state_LQR(:,end+1) = output.value;
    
    % optimal feedback law from DP   
    % INSERT YOU CODE HERE:
    u_DP = ...;
        
    us_DP = [us_DP u_DP];
    cost_DP = [cost_DP cost_DP(end)+u_DP.'*R*u_DP+state_DP(:,end).'*Q*state_DP(:,end)];
    
    % apply DP control to nonlinear system
    input.x = state_DP(:,end);
    input.u = u_DP;
    output = RK4_integrator( @ode, input );
    state_DP(:,end+1) = output.value;
    
    % next time step and visualize result
    iter = iter+1;
    time(end+1) = iter*Ts;
    
    % plot results
    figure(1); clf;
    subplot(221); hold on;
    plot(time,state_LQR(1,:), 'r'); 
    plot(time,state_DP(1,:), 'b'); 
    xlabel('time(s)');
    ylabel('State \phi');
    legend('LQR', 'DP', 'Location', 'northeast')
    
    subplot(223); hold on;
    plot(time,state_LQR(2,:),'r'); 
    plot(time,state_DP(2,:),'b'); 
    xlabel('time(s)');
    ylabel('State \omega');    
    
    subplot(222); hold on;
    stairs(time(1:end-1),us_LQR,'r');
    stairs(time(1:end-1),us_DP,'b');
    xlabel('time(s)');
    ylabel('Control \tau');
    legend('LQR', 'DP', 'Location', 'southeast')
    
    subplot(224); hold on;
    plot(time,cost_LQR,'--ro'); 
    plot(time,cost_DP,'--bo');
    xlabel('time(s)');
    ylabel('Closed-loop cost');
    legend('LQR', 'DP', 'Location', 'southeast')
    
    pause(1e-5);
end

