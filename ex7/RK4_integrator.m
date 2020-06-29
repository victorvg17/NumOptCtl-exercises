function [ output ] = RK4_integrator(ode_fun, input)
% integrator based on the RK4 scheme
% sensitivity generation based on the imaginary trick

% ode_fun is the function handle ode of the system to be integrated, 
% with arguments ode_fun(x, u)
% input is a struct with the following fields:
%  input.x                      initial value
%  input.u                      control (constant)
%  input.Ts                     total integration time
%  input.nSteps                 number of rk4 steps done
%  input.sens                   if true or nonexistent, sensitivies are
%                               computed and returned

% output is the struct output with the following fields:
%  output.value                 output state of integration
%  output.sensX (conditional)   sensitivity of output state wrt to initial
%                               state (partial derivative dxEnd/dx0)
%  output.sensU (conditional)   sensitivity of output state wrt to control
%                               (partial derivative dxEnd/du)

    x0 = input.x;
    u0 = input.u;
    Ts = input.Ts;
    nSteps = input.nSteps;
    
    nx = length(x0);            % state dimension
    nu = length(u0);            % control dimension
    h = Ts / nSteps;            % distance of one rk4 step
    STEP = 1e-100;              % perturbation size for imaginary trick
    
    compute_sensitivities = ~isfield(input,'sens') || input.sens;
    
    xEnd = x0;
    A = eye(nx);
    B = zeros(nx,nu);
    for i = 1:nSteps
        x0 = xEnd;
        xEnd = rk4_step(ode_fun,x0,u0,h);
        if compute_sensitivities
            sensX = zeros(nx,nx); sensU = zeros(nx,nu);
            % imaginary trick for state sensitivity
            for j = 1:nx
                xTemp1 = x0;
                xTemp1(j) = xTemp1(j) + STEP * 1i;
                xTemp1 = rk4_step(ode_fun,xTemp1,u0,h);
                sensX(:,j) = imag(xTemp1)./STEP;
            end
            % imaginary trick for control sensitivity
            for j = 1:nu
                uTemp1 = u0;
                uTemp1(j) = uTemp1(j) + STEP * 1i;
                xTemp1 = rk4_step(ode_fun,x0,uTemp1,h);
                sensU(:,j) = imag(xTemp1)./STEP;
            end
            % propagate sensitivities
            A = sensX*A;
            B = sensX*B + sensU;
        end
    end
    output.value = xEnd;
    if compute_sensitivities
        output.sensX = A;
        output.sensU = B;
    end
end


function x_next = rk4_step(ode_fun,x,u,h)
% one step of the rk4 method
    k1 = ode_fun(x,u);
    k2 = ode_fun(x+h/2.*k1,u);
    k3 = ode_fun(x+h/2.*k2,u);
    k4 = ode_fun(x+h.*k3,u);
    x_next = x + h/6.*(k1+2*k2+2*k3+k4);
end