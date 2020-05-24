function [ P_t ] = ballistic_dynamics_RK4( v_start )
    T = 0.5;
    M = 100;
    DT = T/M;

    X0 = [0, 0, v_start(1), v_start(2), 10, 0, v_start(3), v_start(4)];

    % RK4 integrator
    k1 = ode(X0);
    k2 = ode(X0 + DT*k1'/2);
    k3 = ode(X0 + DT*k2'/2);
    k4 = ode(X0 + DT*k3');
    Xf = X0 + DT/6*(k1' + 2*k2' + 2*k3' + k4');

    P_t = [Xf(1), Xf(2), Xf(5), Xf(6)];
end
