close all
clear variables
clc

HESSIAN_APPROXIMATION = 'exact';
% HESSIAN_APPROXIMATION = 'identity'; rho = 600;
% HESSIAN_APPROXIMATION = 'identity'; rho = 100;
import casadi.*

nv = 2;
x = MX.sym('x',nv);

% Insert your code here (define objective function and its gradient and Hessian
f_expr = .5 * (x(1)-1)^2 + .5 * (10 * (x(2) - x(1)^2))^2 + .5 * x(1)^2;
f = Function('f', {x}, {f_expr});
J = Function('J', {x}, {jacobian(f_expr, x) });
H = Function('H', {x}, {hessian(f_expr, x) });

% Insert your code here (define constraints and its Jacobian and Hessian)
g_expr = x(1) + (1-x(2))^2; 
g = Function('g', {x}, {g_expr});
g_jac = Function('g_jac', {x}, {jacobian(g_expr, x)});
g_hess = Function('g_hess', {x}, {hessian(g_expr, x)});

% SQP solver
max_it = 100;
iter = zeros(nv+1,max_it);
iter(:,1) = [1 -1 1]'; % Initial guess

for i=2:max_it    
    
    x_k = iter(1:2,i-1);
    lambda_k = iter(3,i-1);
    switch HESSIAN_APPROXIMATION
        case 'exact'
        % Insert your code here (exact Hessian)
            Bk = full(H(x_k)) + lambda_k * full(g_hess(x_k)); 
        case 'identity'
        % Insert your code here (scaled identity approximation)
            Bk = rho * eye(2);
    end

    % Build and solve the KKT system
    grad_g = full(g_jac(x_k))'; 
    
    KKT = [Bk, grad_g;
         grad_g', 0 ];
     
    rk = [full(J(x_k))' + lambda_k * grad_g; full(g(x_k))];
    
    iter(:,i) = iter(:,i-1) - KKT \ rk;

end

[X,Y] = meshgrid(-1.5:.05:1.5, -1.5:.05:1.5);
Z = log(1 + 1/2*(X -1).^2 + 1/2*(10*(Y -X.^2)).^2 + 1/2*Y.^2);

y_g = linspace(-0.25,1.5,20);
x_g = -(1 - y_g).^2;

figure(1); clf;
subplot(1,2,1);
surf(X,Y,Z)
hold on;
plot(iter(1,:),iter(2,:), 'ko-')
plot(x_g,y_g,'r');
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
xlabel('x_1')
ylabel('x_2')

subplot(1,2,2); hold on;
plot(iter(1,:),iter(2,:),'ko-')
plot(x_g,y_g,'r');
contour(X,Y,Z);
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
xlabel('x_1')
ylabel('x_2')
legend('solution trajectory','g(x) = 0', 'Location', 'southwest')

figure()
plot(iter(1:2,:)')
grid on
xlabel('iterations')
ylabel('primal solution')
legend('x_1', 'x_2')
grid on