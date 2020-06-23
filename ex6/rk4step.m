function x_next = rk4step(x, u, dynamics, h)
% one rk4 step
% inputs:
%  x             initial state of integration
%  u             control, kept constant over integration
%  dynamics      function handle of ode, callable as dynamics(x, u)
%  h             time step of integration
% output:
%  x_next        state after one rk4 step

    k1 = dynamics(x, u);
    k2 = dynamics(x+h/2.*k1, u);
    k3 = dynamics(x+h/2.*k2, u);
    k4 = dynamics(x+h.*k3, u);
    x_next = x + h/6.*(k1+2*k2+2*k3+k4);
end


