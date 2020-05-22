function [u0,u] = ur(E)
% u''(r)=F(r,E)
% integration: Numerov algorithm
% return u(r) at r=0

rmax = 10;
h = 0.001;
N = rmax / h;
% use array representation
% get the array of F
r = (0:N-1)*h;
F = F_RE(r, E);

% first two points
uStart = rmax*exp(-rmax);
uNext = (rmax-h)*exp(-(rmax-h));

% call Numerov algorithm
u = zeros(1, N);
u = Numerov(-h, N, 2, F, uStart, uNext, u);
u(1)=2*u(2)-u(3)+h^2*F(2)*u(2);

u0 = u(1);

end