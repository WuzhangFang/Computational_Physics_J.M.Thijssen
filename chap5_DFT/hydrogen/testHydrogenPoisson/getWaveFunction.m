function [u0,u] = getWaveFunction(E,V,h,N)
    % first two points of radial wave function
    rmax=N*h;
    r=(0:N-1)*h;
    uStart = rmax*exp(-rmax);
    uNext = (rmax-h)*exp(-(rmax-h));
    F = F_RE(r,E) + 2*V;
    % call Numerov algorithm
    u = zeros(1, N);
    u = Numerov(-h, N, 2, F, uStart, uNext, u);
    u(1)=2*u(2)-u(3)+h^2*F(2)*u(2);
    u0 = u(1);
end
