function lamJ = back_AD(u,m,param)

N = length(u);
h = param.h;
x0 = param.x0;

% initialize bar quantities
ubar = zeros(N,1);
xbar = zeros(N,1);
xbar(m) = 1;

% forward sweep
X = zeros(N,1);
X(1) = x0;
for i = 1:N-1
    X(i+1) = X(i) + h * ((1 - X(i)) * X(i) + u(i));
end

% backward sweep
for i = N:-1:2
   dfdu = h;
   dfdx = 1 + h * (1 - 2 * X(i));
   ubar(i) = ubar(i) + xbar(i) * dfdu; 
   xbar(i-1) = xbar(i-1) + xbar(i) * dfdx; 
end
% x1 = f(x0,u)
ubar(1) = ubar(1) + xbar(1) * h;

% mth row of jacobian
lamJ = ubar';