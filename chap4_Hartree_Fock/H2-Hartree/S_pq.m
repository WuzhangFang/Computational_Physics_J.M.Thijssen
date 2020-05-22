function y = S_pq(p, Rp, q, Rq)
% calculates the element of overlap matrix S_pq
R = Rp - Rq;
K = exp(-p*q*R^2/(p+q));
y = (pi/(p+q))^1.5*K;

end