function [u0,u,n] = getWaveFunction(E,U,h,N)
    % u0: u(r=0)
    % u: radial wave function
    % n: density
    rmax=N*h;
    uStart = rmax*exp(-2*rmax);
    uNext = (rmax-h)*exp(-2*(rmax-h));
    F = F_RE(E,U,h,N);
    % call Numerov algorithm
    u = zeros(1, N);
    u = Numerov(-h, N, 2, F, uStart, uNext, u);
    u(1)=2*u(2)-u(3)+h^2*F(2)*u(2);
    u0 = u(1);
    n=zeros(1,N);
    for i=1:N
        n(i)=u(i)^2;
    end
    n(1)=0;

    normFactor=1/calInt(n,N,h);
    n(1)=0;
    for i=2:N
        R=(i-1)*h;
        n(i)=-normFactor*n(i)/R;
    end
end
