% solve U''(r)=-u(r)^2/r (5.77)
% use array representation

rmax=10;
h=0.01;
N=rmax/h;

% get the array of F
r = (0:N-1)*h;

% first two points
UStart = 0;
UNext = h;

% normalize wave function u(r)
[~,u] = getWaveFunction(-0.5,0,h,N);
u = normalize(u,r);
%plot(r,u);


n = -u.^2./r; % (5.77)
n(1)=0;
% plot(r,n);

% call Verlet algorithm
U = zeros(1,N);
U = NumInhom(h,0,rmax,UStart,UNext,n,U);

% determin alpha
zScr=1-(rmax^2+2*rmax+1)*exp(-rmax);
alpha = (zScr - U(N))/ (rmax-h);

U = alpha*r + U;

% plotting
Uanaly = -(r+1).*exp(-2*r)+1; % (5.81)
plot(r,U,r,Uanaly);
legend('Numerov','Analytical');
