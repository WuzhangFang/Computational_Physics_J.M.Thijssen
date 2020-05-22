function y = T_pq(p, Rp, q, Rq)
% calculates the element of kinetic matrix T_pq      
R = Rp - Rq;
K = exp(-p*q*R^2/(p+q));
y = (p*q)*pi^1.5*(3-2*p*q*R^2/(p+q))*K /(p+q)^2.5;

end