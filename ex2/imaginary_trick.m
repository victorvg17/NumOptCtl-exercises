F = @(x) x.^2;
t = 10^0;
u = 1;
F_der = F(u + 1i*t*u);
F_der = imag(F_der)/t