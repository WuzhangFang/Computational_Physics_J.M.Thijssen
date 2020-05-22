function u = normalize(u, r)
% r and u are arrays
% y: the nomralizetion factor
f = trapz(r, u.^2);
u = u / sqrt(f);
end