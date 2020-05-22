function [r1, r2, u1, u2] = Numerov(L, E, h, rmax)
% h: integration interval
% r1 and r2 are the positions beyond rmax
% L: angular momentum quantum number

% parameters
alpha=6.12; % scaling of h^2/2m
epsilon=5.9; % meV
C=sqrt(epsilon*alpha/25);
r0=0.6; % starting point

% initialization
n1 = floor((rmax-r0)/h) + 1;
n2 = n1 + floor((2*pi/sqrt(alpha*E)) / h);
r = [0 h (2:n2)*h]+r0;
u = zeros(1, n2);
w = zeros(1, n2);

% first point
u(1) = exp(-C*r0^(-5)); % (2.17)
u_dot = exp(-C*r0^(-5))*5*C*r0^(-6);
w(1) = (1 - (h^2/12)*F_LRE(L,r0,E))*u(1); % (2.13)

% second point
numerator = (2+5*h^2*F_LRE(L,r0,E)/6)*(1-h^2*F_LRE(L,r0-h,E)/12)*u(1)+2*h*u_dot*(1-h^2*F_LRE(L,r0-h,E)/6);
denominator = (1-h^2*F_LRE(L,r0+h,E)/12)*(1-h^2*F_LRE(L,r0-h,E)/6)+(1-h^2*F_LRE(L,r0-h,E)/12)*(1-h^2*F_LRE(L,r0+h,E)/6);
u(2) = numerator / denominator;
w(2) = (1 - (h^2/12)*F_LRE(L,r0+h,E))*u(2); % (2.13)

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