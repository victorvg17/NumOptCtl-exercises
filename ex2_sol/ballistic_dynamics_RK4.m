function [ P_t ] = ballistic_dynamics_RK4( v_start )

    T = 0.5;
    M = 100;
    DT = T/M;

    X0 = [0, 0, v_start(1), v_start(2), 10, 0, v_start(3), v_start(4)];
    Xf = X0;
    % RK4 integrator
    for j=1:M
        % insert your code calling the provided ode() function here
        Xf = rk4step(@(t,x) ode(x), DT, Xf, 0 ) ;
    end

    P_t = [Xf(1), Xf(2), Xf(5), Xf(6)];

end

