function sigmaL = crossSection(L,E,h,rmax)
% return the cross section of L
% energy in meV, radius in rho=3.57 angstrom
[r1,r2,u1,u2] = Numerov(L,E,h,rmax);
alpha=6.12; % scaling of h^2/2m
k = sqrt(alpha*E);
K = (r1*u2)/(r2*u1); %(2.9b)
[j1, n1] = Bessel(L,k*r1);
[j2, n2] = Bessel(L,k*r2);
delta = atan((K*j1-j2)/(K*n1-n2)); %(2.9a)
sigmaL = ((4*pi)/k^2) * (2*L+1)*sin(delta)^2; %(2.8)
end