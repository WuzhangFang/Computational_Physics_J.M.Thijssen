function y = F_LRE(L, r, E)
% energy in meV, radius in rho=3.57 angstrom
alpha=6.12; % scaling of h^2/2m
y = alpha * (V(r) + L*(L+1)/(alpha*r^2) - E);

end