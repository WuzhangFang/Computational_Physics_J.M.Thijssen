%test harmonic potential
L=0;
E=3;
h = 0.01;
rmax=1;

A = exp(h^2/2); % very close to 1
[r1, r2, u1, u2] = Numerov(L,E,h,rmax);

fprintf('Analytical u(r1): %f\n', A*r1*exp(-r1^2/2));
fprintf('Numerov result u(r1): %f\n', u1);
