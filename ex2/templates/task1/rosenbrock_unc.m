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

%% 1b) implement gradient and hessian
clc; clear; close all;
% grad = @(x, y) [2*(x - 1);
%                 40*(x - 1) + 200*(y - x.^2)];
% hess = @(x, y) [[2, 40*(1 - 10*x)];
%                 [0, 200]];

grad = @(x, y) [2*(x-1) - 400*x*y + 400*x.^3;
                200*(y-x.^2)];
hess = @(x, y) [[2 - 400*y + 1200*x.^2, -400*x];
                [-400*x, 200]];

            
N = 1000;
w0 = [1; 1.1];
%exact-Newton type
w_ex_nt = w0;
for k = 2:N
   w_ex_nt(:, k) = newton_step(hess, grad, w_ex_nt(:, k-1));
end

figure(1); clf; hold on;
plot(w_ex_nt(1, :), w_ex_nt(2, :), '-+', 'DisplayName', 'exact Newton');
legend('Location', 'east outside');
xlabel('x'); ylabel('y');

%% 1c)approximate Hessian
rho = [100, 500];
hess_aprx = @(x, y) [[rho(2), 0];
                     [0, rho(2)]];
w_aprx_nt = w0;
for k = 2:N
   w_aprx_nt(:, k) = newton_step(hess_aprx, grad, w_aprx_nt(:, k-1));
end
plot(w_aprx_nt(1, :), w_aprx_nt(2, :), '-+', 'DisplayName', 'apprx Newton');

%% 1d) using casadi
w_ca = MX.sym('w_ca', 2, 1);
expr = (1 - w_ca(1)).^2 + 100*(w_ca(2) - w_ca(1).^2).^2;
j_expr = jacobian(expr, w_ca);
J = Function('J', {w_ca}, {j_expr});
h_expr = hessian(expr, w_ca);
H = Function('H', {w_ca}, {h_expr});

w_ca = w0;
for k = 2:N
    hess = H(w_ca(:, k-1));
    jak = J(w_ca(:, k-1));
    w_ca(:, k) = full(w_ca(:, k-1) - inv(hess)*jak');
end
plot(w_ca(1, :), w_ca(2, :), '-+', 'DisplayName', 'casadi newton');