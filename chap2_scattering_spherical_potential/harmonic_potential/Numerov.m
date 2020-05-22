function [r1, r2, u1, u2] = Numerov(L, E, h, rmax)
% h: integration interval
% r1 and r2 are the positions beyond rmax
% L: angular momentum quantum number

n1 = floor(rmax / h) + 1;
n2 = n1 + floor((2*pi/sqrt(E)) / h);
r = (1:n2)*h;

% initialization
u = zeros(1, n2);
u(1) = 0; % u(0) = 0
u(2) = h^(L+1); % u(h)=h^(L+1)

w = zeros(1, n2);
w(1) = 0; % w(0) = 0
w(2) = (1 - (h^2/12)*F_LRE(L,r(2),E))*u(2); % from (2.13)

% integration to r1 and r2
for i=3:n2
w(i) = 2*w(i-1) - w(i-2) + (h^2)*F_LRE(L,r(i-1),E)*u(i-1);
u(i) = (1 - (h^2/12)*F_LRE(L,r(i),E))^(-1)*w(i);
end


% return values
r1 = r(n1);
r2 = r(n2);
u1 = u(n1);
u2 = u(n2);

end