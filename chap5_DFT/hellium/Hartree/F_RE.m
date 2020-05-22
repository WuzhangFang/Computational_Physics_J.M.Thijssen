function F = F_RE(E, U, h, N)
% u''(r)=F(r)
% r: array
% For U=0, agrees with Fortran code
F = zeros(1,N);

for i=2:N
    R = (i-1)*h;
    F(i) = -2 * (E + (2 - U(i))/R);
end

end