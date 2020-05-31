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
import casadi.*
 
% Create NLP: Solve the Rosenbrock problem:
%     minimize    (1-x)^2 + 100*z^2
%     subject to  z + x^2 - y == 0
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
v = [x;y;z];
f = (1-x)^2 + 100*z^2;
g = z - (y - x^2);
nlp = struct('x', v, 'f', f', 'g', g);

% Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);

% Solve the NLP
res = solver('x0' , [2.5 3.0 0.75],... % solution guess
             'lbx', -inf,...           % lower bound on x
             'ubx',  inf,...           % upper bound on x
             'lbg',    0,...           % lower bound on g
             'ubg',    0);             % upper bound on g
 
% Print the solution
f_opt = full(res.f)          
x_opt = full(res.x)          
lam_x_opt = full(res.lam_x) 
lam_g_opt = full(res.lam_g)  
%% 2b)implement and solve optimizatio problem
clc; clear; close all;
v_bar = 15;
% alpha = -1:0.2:1;
alpha = -0.4;
v = MX.sym('v', 2, 1);
p = ballistic_dynamics_RK4(v);
f = -p(1);
g = [-p(2);
    -alpha*(p(1)-10)-p(2);
    v(1)^2 + v(2)^2];

nlp = struct('x', v, 'f', f', 'g', g);
solver = nlpsol('solver', 'ipopt', nlp);

% Solve the nlp
res = solver('x0', [0, 0], ... % solution guess
             'lbx', -inf,  ... % lower bound on x
             'ubx', inf,   ... % upper bound on x
             'lbg', [-inf; -inf; -inf], ... % lower bound on g
             'ubg', [0; 0; v_bar^2]);

% Print the solution
f_opt_ball = full(res.f);
x_opt_ball = full(res.x);
lam_x_opt_ball = full(res.lam_x);
lam_g_opt_ball = full(res.lam_g);