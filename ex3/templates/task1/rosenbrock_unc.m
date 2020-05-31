%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
% Load CasADi
clear variables
close all
import casadi.*
 
% Create NLP: Solve the Rosenbrock problem:
%     minimize    (x-1)^2 + 100*(y - x^2)^2
%

x = SX.sym('x');
y = SX.sym('y');
v = [x;y];
f = unc_rosenbrock_fun(x,y);
g = [];
nlp = struct('x', v, 'f', f', 'g', g);

% Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);

% % Solve the NLP
res = solver('x0' , [2.5 3.0],...      % solution guess
             'lbx', -inf,...           % lower bound on x
             'ubx',  inf);             % upper bound on g
              
 
% Print the solution
f_opt = full(res.f)          
x_opt = full(res.x)         
lam_x_opt = full(res.lam_x)  
lam_g_opt = full(res.lam_g)  

% Plot result
figure()
n_points = 50;
range = 2;

x_v = linspace(x_opt(1) - 2*range,  x_opt(1) + range,n_points);
y_v = linspace(x_opt(2) - 2*range,  x_opt(2) + 2*range,n_points);

[X, Y] = meshgrid(x_v, y_v);

surf(X, Y, unc_rosenbrock_fun(X,Y) )
xlabel('x')
ylabel('y')
zlabel('z')
%% 1c) implement f, g, their jacobians and hessians
clc; clear; close all;
f_obj = @(x) 1/2*(x(1)-1)^2 + ...
             1/2 * (10 * (x(2)-x(1)^2) )^2 + 1/2*x(1)^2;
g_constr = @(x) x(1) + (1-x(2))^2;

f_jac = @(x) 200*[x(1)^3-x(1)*x(2)+x(1)/200;
                  1/2*(x(2)-x(1)^2)];
g_jac = @(x) [1;
              -2*(1-x(2))];

f_hess = @(x) 200*[3*x(1)^2-x(2)+1/200, -x(1);
                   -x(1), 1/2];
g_hess = @(x) [0, 0;
               0, 2];

x = MX.sym('x', 2, 1);
f = Function('f', {x}, {f_obj(x)});
f = returntypes('full', f);

g = Function('g', {x}, {g_constr(x)});
g = returntypes('full', g);

delta_f = Function('delta_f', {x}, {f_jac(x)});
delta_f = returntypes('full', delta_f);

delta_g = Function('delta_g', {x}, {g_jac(x)});
delta_g = returntypes('full', delta_g);

hess_f = Function('hess_f', {x}, {f_hess(x)});
hess_f = returntypes('full', hess_f);

hess_g = Function('hess_g', {x}, {g_hess(x)});
hess_g = returntypes('full', hess_g);

%% 1d)implement Newton type method
%% i) constant Hessian apprx
rho = 0:100:600;
B = 1*eye(2);
w0 = [1; -1; 1];
w = w0;

for i = 1:length(rho)
    Bi = rho(i)*eye(2);
    for k = 2:100
        xk = w(:, k-1);
        Mk = [Bi, delta_g([xk(1); xk(2)]);
              delta_g([xk(1); xk(2)])', 0];
        delta_fk = delta_f([xk(1), xk(2)]);
        val_gk = g([xk(1), xk(2)]);
        Ak = [delta_fk; val_gk];
        w(:, k) = w(:, k-1) - pinv(Mk)*Ak;
    end
    W_mother(:, :, i) = w;
end

figure(1); clf; hold on;
wi = W_mother(:, :, 1);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:0');
wi = W_mother(:, :, 2);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:100');
wi = W_mother(:, :, 3);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:200');
wi = W_mother(:, :, 4);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:300');
wi = W_mother(:, :, 5);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:400');
wi = W_mother(:, :, 6);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:500');
wi = W_mother(:, :, 7);
plot(wi(1,:), wi(2,:), '-+', 'DisplayName', 'rho:600');
legend('Location', 'east outside');
xlabel('x');
ylabel('y');
hold off;

%% ii) Lagrangian Hessian apprx
figure(2); hold on;
w_exact = w0;
for k = 2:100
    xk = w_exact(:, k-1);
    Bk = hess_f([xk(1); xk(2)]) + hess_g([xk(1); xk(2)]);
    Mk = [Bk, delta_g([xk(1); xk(2)]);
          delta_g([xk(1); xk(2)])', 0];
    delta_fk = delta_f([xk(1), xk(2)]);
    val_gk = g([xk(1), xk(2)]);
    Ak = [delta_fk; val_gk];
    w_exact(:, k) = w_exact(:, k-1) - pinv(Mk)*Ak;
end
plot(w(1,:), w(2,:), '-+', 'DisplayName', 'exact Hessian');
legend('Location', 'east outside');
xlabel('x');
ylabel('y');