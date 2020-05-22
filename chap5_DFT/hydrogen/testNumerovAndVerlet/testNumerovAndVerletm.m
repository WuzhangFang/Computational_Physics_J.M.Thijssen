% solve the radial equation of hydrogen (5.72) 
% 1. Numberov algorithm
% 2. Verlet algorithm

rmax = 10;
h = 0.001;
N = rmax / h;
% initialization
% F = zeros(1, N);
E = -0.5; % bound state of hydrogen
% use array representation
% get the array of F
r = (0:N-1)*h;
F = F_RE(r, E);
% for i=2:N
%     r = (i-1) * h;
%     F(i) = F_RE(r, E);
% end

% first two points for Numerov algorithm
uStart = rmax*exp(-rmax);
uNext = (rmax-h)*exp(-(rmax-h));

% call Numerov algorithm
u = zeros(1, N);
u = Numerov(-h, N, 2, F, uStart, uNext, u);
% u(r)=rR(r) at r=0 R(r) is diverged, assume w(r)=u(r)
u(1)=2*u(2)-u(3)+h^2*F(2)*u(2);

% call Verlet algorithm
v = zeros(1, N);
v = Verlet(-h, N, 1, F, uStart, uNext, v);
%v(1)=2*v(2)-v(3)+h^2*F(2)*v(2);

% plotting along with analytical formula u(r)=r*exp(-r)
plot(r,u,'b--',r,v,'r:');
legend('Numerov','Verlet');
