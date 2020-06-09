clear; clc; close all;
f = @(x, y) (1-x).^2 + 100*(y - x.^2).^2;
gradf = @(x, y) [2*(x-1) - 400*x*y + 400*x.^3;
                200*(y-x.^2)];
hessf = @(x, y) [2 - 400*y + 1200*x.^2, -400*x;
                -400*x, 200];

grad = @(w) gradf(w(1), w(2));
w0 = [1; 1.1];
numit = 1000;

figure(1); clf; hold on;
x = linspace(0.8, 1.2, 100);
[X,Y] = meshgrid(x);
contour(X,Y, f(X,Y), 100, 'DisplayName', 'f(w)');

scatter(w0(1), w0(2), 'c*', 'DisplayName', 'w_0')
axis([min(x), max(x), min(x), max(x)])
xlabel('x')
ylabel('y')
legend()

%% gradient descent , rho = 100
W_gd100 = newton_type(w0, grad, @(w) 100 * eye(2), numit);
plot(W_gd100(1,:), W_gd100(2,:), 'rx-', 'DisplayName', 'GD, \rho=100' )

%% gradien desc, rho = 500
W_gd500 = newton_type(w0, grad, @(w) 500* eye(2), numit);
plot( W_gd500(1,:), W_gd500(2,:), 'bo-', 'DisplayName', 'GD, \rho=500' )

%% exact hessian
W_eh = newton_type(w0, grad, @(w) hessf(w(1), w(2)), numit);
plot( W_eh(1,:), W_eh(2,:), 'k^-', 'DisplayName', 'exact Hessian' )

%% plot of convergence speed
figure(2); clf; hold on;
xopt = [1;1];                               % true minimizer
plot(0:numit, log10(max(abs(W_gd500 - xopt))))
plot(0:numit, log10(max(abs(W_eh - xopt))))
legend('GD, \rho=500', 'exact hessian')
xlabel('iteration k')
ylabel('log ||w_k - w^* ||_\infty')

%% Casadi
import casadi.*
w = MX.sym('w', 2);
f_expr = f(w(1), w(2));
f_cas = Function('f_cas', {w}, {f_expr});
grad_cas = Function('grad_cas', {w}, { jacobian(f_expr, w)'});
grad_cas = returntypes('full', grad_cas);
hess_cas = Function('hess_cas', {w}, { hessian(f_expr, w) });
hess_cas = returntypes('full', hess_cas);

% exact hessian
W_cas = newton_type(w0, grad_cas, hess_cas, numit);

figure(1)
plot(W_cas(1,:), W_cas(2,:), 'gv-','DisplayName', 'exact Hessian, casadi');


