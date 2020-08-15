function Jp = forw_AD(u,m,param)

N = length(u);
h = param.h;
x0 = param.x0;

udot = zeros(N,1);
udot(m) = 1;
xdot = zeros(N,1);

xk = x0;
xdot(1) = h * udot(1);

for i = 2:N
    xk = xk + h * ((1 - xk) * xk + u(i-1)); 
    dfdu = h;
    dfdx = 1 + h * (1 - 2*xk);
    xdot(i) = dfdu * udot(i) + dfdx * xdot(i-1);
end

% m-th column of jacobian
Jp = xdot;