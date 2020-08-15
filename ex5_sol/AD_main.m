close all; clear; clc;
import casadi.*

%%
x0 = 0.5;
N = 500;
h = .01;

param.x0 = x0;
param.h = h;

x = SX.sym('x',1);
u = SX.sym('u',1);
U = SX.sym('U',N);

f = Function('f',{x,u}, {x + h * ((1 - x) * x + u)});

Phi_expr = x0;
for i = 1:N
    Phi_expr = [Phi_expr; f(Phi_expr(i), U(i))];
end
Phi_expr = Phi_expr(2:end);

Phi = Function('Phi', {U}, {Phi_expr});
J = Function('J', {U}, {jacobian(Phi_expr, U)} );

%% tests

utest = rand(N,1);
Jref = full(J(utest));
m = 1;

fprintf('Mismatch between hand-coded forward AD and CasADi (column %d):\n',m)
disp(max(max(abs(Jref(:,m) - forw_AD(utest, m, param)))))

m = N;
fprintf('Mismatch between hand-coded backward AD and CasADi (row %d):\n',m)
disp(max(max(abs(Jref(m,:) - back_AD(utest, m, param)))))

fprintf('Mismatch between hand-coded forward AD and CasADi (full jacobian)\n')
disp(max(max(Jref(N-m+1:end,:) - J_FAD(utest, m, param))))

fprintf('Mismatch between hand-coded backward AD and CasADi (full jacobian)\n')
disp(max(max(Jref(N-m+1:end,:) - J_BAD(utest, m, param))))



%% timing

t_fw = zeros(N,1);
t_bw = zeros(N,1);

for i = 1:N
    tic
    J_FAD(utest, i, param);
    t_fw(i) = toc;
    
    tic
    J_BAD(utest, i, param);
    t_bw(i) = toc;
end

figure(1); clf; hold on;
plot(t_fw);
plot(t_bw);
legend('FAD', 'BAD')
xlabel('rows computed')
ylabel('t in s')