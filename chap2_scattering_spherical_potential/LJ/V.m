% Lennard-Jones potential (2.14)
function v = V(r)
% energy in meV, radius in rho=3.57 angstrom
rmax = 5; % need to adjust
epsilon=5.9; % meV
if r <= rmax
    v = epsilon * (1./r.^12-2./r.^6);
else
    v =0;
end