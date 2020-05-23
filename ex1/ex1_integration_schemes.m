Ts = 0.5;
tspan = [0*Ts 19*Ts];
x0 = [1; 0];
h_eul = 0.125;
h_rk = 0.5;
[t,x] = ode45(@derivative_fun,tspan, x0);

[t_eul, x_eul] = euler_integrator(tspan, x0, h_eul);
[t_rk, x_rk] = rk_integrator(tspan, x0, h_rk);

% plot results
figure(1)
plot(t, x(:, 1), '-o', t, x(:, 2), '-o');

hold on 
plot(t_eul, x_eul(:, 1), '-+', t_eul, x_eul(:, 2), '-+');
plot(t_rk, x_rk(:, 1), '-s', t_rk, x_rk(:, 2), '-s');
title('damped oscillator solution using ode45 solver');
xlabel('Time');
ylabel('Solution X');
legend('ode1', 'ode2', 'eul1', 'eul2', 'rk1', 'rk2');

hold off

%% rk integrator for lorenz
param.rho = 28;
param.beta = 8/3;
param.sigma = 10;
x0 = [1; 0; 0];
t_span = [0 100];
h = 0.01;

[t_rk, x_rk] = rk_integrator_lorenz(t_span, x0, h, param);

figure(2)
plot3(x_rk(:, 1), x_rk(:, 2), x_rk(:, 3), '-o');
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('lorenz attractor');