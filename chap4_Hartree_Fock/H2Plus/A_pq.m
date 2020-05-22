function y = A_pq(p, Ra, q, Rb, Rc)
% calculates the element of coulomb matrix A_pq
% p, Ra: basis function p at Ra
% q, Rb: basis function q at Rb
% Rc: Coulomb potential 1/(r-Rc)

K = exp(-p*q*(Ra-Rb)^2/(p+q));
Rp = (p*Ra + q*Rb) / (p+q);

y = (-2*pi*K*F0((p+q)*(Rp-Rc)^2))/(p+q);
end