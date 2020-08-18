% Implementation of an interior point method 
clear variables
close all
clc

import casadi.*
% Problem definition
nv = 2;
ne = 1;
ni = 1;

x = MX.sym('x',nv);

x_test = [2,3];
numtol = 1e-12;
% YOUR CODE HERE:
% Inequality contstraints
% H = Function(...);
h_test = 9.0;
H = Function('H',{x}, {x(1)^2 + x(2)^2 - 4});
if abs(full(H(x_test)) - h_test) >= numtol
    error('error in H')
end

% Equality contstraints
% G = Function(...);
G_test = -8.090702573174319;
G = Function('G',{x}, {sin(x(1)) - x(2)^2});
if abs(full(G(x_test)) - G_test) >= numtol
    error('error in G')
end

% Objective
% F = Function(...);
F_test = 5;
F = Function('F',{x}, {(x(1) - 4)^2 + (x(2) - 4)^2 });
if abs(full(F(x_test)) - F_test) >= numtol
    error('error in F')
end

% Jacobian of g
% Jg = Function(...);
Jg_test = [-0.4161,   -6.0000];
Jg = Function('Jg',{x},{ jacobian(G(x), x) });
if min(abs(full(Jg(x_test)) - Jg_test)) >= numtol
    error('error in Jg')
end


% Jacobian of h
% Jh = Function(...);
Jh_test = [4,   6];
Jh = Function('Jh',{x},{ jacobian(H(x), x) });
if min(abs(full(Jh(x_test)) - Jh_test)) >= numtol
    error('error in Jh')
end


% Jacobian of f
% Jf = Function(...);
Jf_test = [4,   -2];
Jf = Function('Jf',{x},{ jacobian(F(x), x) });
if min(abs(full(Jf(x_test)) - Jf_test)) >= numtol
    error('error in Jf')
end


% Hessian of g
% Hg = Function(...);
Hg_test =    [  -0.9093,  0        ; ...
                0,        -2.0000 ];
Hg = Function('Hg',{x},{hessian(G(x), x)});
if min(min(abs(full(Hg(x_test)) - Hg_test))) >= numtol
    error('error in Hg')
end
            
% Hessian of h
% Hh = Function(...);
Hh_test =   [   2, 0; ...
                0  2];
Hh = Function('Hh',{x},{hessian(H(x), x)});
if min(min(abs(full(Hh(x_test)) - Hh_test))) >= numtol
    error('error in Hh')
end
% Hessian of f
% Hf = Function(...);
Hf_test =   [   2, 0; ...
                0  2];
Hf = Function('Hf',{x},{hessian(F(x), x)});
if min(min(abs(full(Hf(x_test)) - Hf_test))) >= numtol
    error('error in Hf')
end

% Interior point solver
max_it = 100;
% xk = [-2; 4];
xk = [-2; -4];
lk = 10*ones(ne,1);
vk = 10*ones(ni,1);
sk = 10*ones(ni,1);
iter = zeros(nv + ne + ni + ni,max_it);
iter(:,1) = [xk;lk;vk;sk];
tau = 2;
k_b = 1/3;
th_1 = 1.0e-8;
th_2 = 1.0e-8;
for i = 2:max_it
    % Build KKT system
    Hf_e    = Hf(xk);
    Hg_e    = Hg(xk);
    Hh_e    = Hh(xk);
    Hl      = Hf_e + Hg_e*lk + Hh_e*vk;
    
    Jg_e    = Jg(xk);
    Jh_e    = Jh(xk);
    Jf_e    = Jf(xk);
    
    g_e     = G(xk);
    h_e     = H(xk);
    
    % YOUR CODE HERE:
    % Buiild the KKT system
    M = full([Hl    Jg_e'  Jh_e'        zeros(nv, ni);   
              Jg_e    0     0           0;   
              Jh_e    0     0           eye(ni);
              0 0     0     diag(sk)     diag(vk)]);
    
    rhs = - full([ Jf_e' + Jg_e' * lk + Jh_e' * vk ;
                    g_e;
                    h_e + sk;
                    vk*sk - tau]);
    
    % Termination condition
    if norm(rhs) < th_1
        if tau < th_2
            display('Solution found!')
            break;
        else
            tau = tau*k_b;
        end
    end
     
    % YOUR CODE HERE:
    % Compute Newton step
    sol = M \ rhs;
    
    % line-serach
    max_ls = 100;
    x_step  = sol(1:nv);
    l_step  = sol(nv+1:nv+ne);
    v_step  = sol(nv+ne+1:nv+ne+ni);
    s_step  = sol(nv+ne+ni+1:end);
    alpha = 1;
    k_ls = 0.9;
    min_step = 1.0e-8;
    for j=1:max_ls
        
        % YOUR CODE HERE: 
        % Compute trial step
        v_t = vk + alpha * v_step;
        s_t = sk + alpha * s_step;
        if all(v_t >= 0) && all(s_t >= 0)
            break;
        end
        
        % YOUR CODE HERE:
        % Decrease alpha
        alpha = k_ls * alpha;
        
        % Terminiation condition
        if norm(alpha*[ v_step;s_step]) < min_step
            error('Line search failed! Could not find dual feasible step.')
        end
    end
    
    % actual step
    xk  = xk + alpha*x_step;
    lk  = lk + alpha*l_step;
    vk  = vk + alpha*v_step;
    sk  = sk + alpha*s_step;
    % save for later processing
    iter(:,i) = [xk;lk;vk;sk];

    % Print some results
    if mod(i,20) == 1           % every now and then reprint header
        fprintf([repmat('-',1,49),'\n']);
        fprintf('it \t tau \t\t ||rhs|| \t alpha\n');
        fprintf([repmat('-',1,49),'\n']);
    end
    fprintf('%d \t %e \t %e \t %e\n',i, tau, norm(rhs), alpha)

end
iter = iter(:,1:i-1);
plot(iter')
grid on
xlabel('iterations')
ylabel('solution')
legend('x1', 'x2' ,'\lambda' ,'\nu', 's')

% Plot feasible set, and iterations
figure()
pl = ezplot('sin(x) - y^2');
set(pl,'Color','red');
hold all
pl= ezplot('x^2 + y^2 - 4');
set(pl,'Color','blue');
ezcontour('(x- 4)^2 + (y - 4)^2')
plot(iter(1,:),iter(2,:),'--','Color','black')
plot(iter(1,:),iter(2,:),'o','Color','black')
title('Iterations in the primal space')
grid on

%% check SOSC
tol = 1e-6;                 % tolerance for when sth is zero

strict_complementarity = true;
% check for activness of h
if abs(full(h_e)) <= tol               % h is active
    Jg_tilde = full([Jg_e; Jh_e]);     % extended equality constr jacobian
    if vk >= tol
        disp('h is strictly active')
    else
        disp('h is weakly active')
        strict_complementarity = false;
    end
else
    disp('h is inactive');
    Jg_tilde  = full(Jg_e);
end

% compute reduced Hessian
Z = null(Jg_tilde);         % nullspace
redH = Z' * full(Hl) * Z;         % reduced Hessian
eigs = eig(redH);
mineig = min(eigs);

if ~strict_complementarity
    disp('strict complentarity does not hold.')
    disp('The conditions for the the theorem second order optimality conditions are not fulfilled.')
elseif isempty(Z) ||  mineig > tol
    disp('redH > 0. SOSC (and SONC) is fullfilled');
    disp('The solution is a local minimizer.');
elseif mineig >= -tol
    disp('redH >= 0. SONC is fullfilled')
    disp('The solution might be a local minimizer.');
else
    disp('redH not PSD. Neither SONC nor SOSC hold.')
    disp('The solution is not a local minimizer.');
end
