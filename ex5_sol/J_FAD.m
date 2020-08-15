function J = J_FAD(u,m,param)
N = length(u);
h = param.h;
x0 = param.x0;

J = zeros(m,N);

% forward sweep
X = zeros(N,1);
X(1) = x0;
for i = 1:N-1
    X(i+1) = X(i) + h * ((1 - X(i)) * X(i) + u(i));
end

% forward AD
for n = 1:N
    udot = zeros(N,1);
    udot(n) = 1;
    xdot = zeros(N,1);

    xdot(1) = h * udot(1);
    
    for i = 2:N
        dfdu = h;
        dfdx = 1 + h * (1 - 2*X(i));
        xdot(i) = dfdu * udot(i) + dfdx * xdot(i-1);
    end
    J(:,n) = xdot(N-m+1:N);
end