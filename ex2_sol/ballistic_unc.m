clc; clear; close all;
import casadi.*


% decision variables
alpha1 = SX.sym('alpha1',1);
alpha2 = SX.sym('alpha2',1);
w = [alpha1; alpha2];

vbar = 15;
v = @(alpha)  vbar * [cos(alpha); sin(alpha)];

% expression for objective
P =  ballistic_dynamics_RK4([ v(alpha1) ; v(alpha2)])';
p1 = P(1:2);
p2 = P(3:4);
dp = p1 - p2;
f_expr = dp' * dp;

% gradient and hessian
f = Function('f', {w}, {f_expr});
grad_f = Function('grad_f', {w}, {jacobian(f_expr,w)'});
grad_f = returntypes('full',grad_f);
hess_f = Function('hess_f', {w}, {hessian(f_expr,w)});
hess_f = returntypes('full',hess_f);

%% solve
numit = 20;
w0 = [pi/4; 3 * pi/4];

W = newton_type(w0, grad_f, hess_f, numit);
wopt = W(:,end);
disp(['objective value at sol: ', num2str(full(f(wopt)))]);

figure(1); clf; hold on;
plot(0:numit, W(1,:))
plot(0:numit, W(2,:))
xlabel('iteration k')
legend('\alpha_1', '\alpha_2')

% simulate optimal solution

v1opt = v(wopt(1));
v2opt = v(wopt(2));

T = 0.5;
M = 100;
DT = T/M;
X0 = [0, 0, v1opt', 10, 0, v2opt']';
X = X0;

% RK4 integrator
for j=1:M
    % insert your code calling the provided ode() function here
    X(:,j+1) = rk4step(@(t,x) ode(x), DT, X(:,j), 0 ) ;
end

figure(2); clf; hold on;
plot(X(1,:), X(2,:))
plot(X(5,:), X(6,:))
