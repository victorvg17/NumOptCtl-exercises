function LQR_design(Q,R,Ts,nSteps)
% creates LQR feedback matrix for linearized discrete time system 
%  and saves it in lqr.mat
% inputs:
%   Q           weighting matrix for states
%   R           weighting matrix for controls
%   Ts          time discretization constant
%   nSteps      number of integration steps during Ts

% linearize at x_lin, u_lin
x_lin = [0; 0];
u_lin = 0;
input.Ts = Ts;
input.nSteps = nSteps;
input.x = x_lin;
input.u = u_lin;
% ode of the systems
% ode = @ (x, u)[x(2); sin(x(1) + u)];
output = RK4_integrator( @ode, input );

% matrices of linearized system
% x_k+1 = A * x_k + B * u_k
A = output.sensX;
B = output.sensU;

% design LQR controller
% INSERT YOUR CODE HERE:
[K,P] = dlqr(A,B,Q,R);

save lqr.mat A B Q R K P
end