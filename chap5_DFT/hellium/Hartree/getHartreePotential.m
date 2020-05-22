function U = getHartreePotential(n,h,N)
    % solve U''(r)=n
    % n: density
    % U: Hartree potential

    rmax = N*h;
    % call Numerov algorithm to solve U''(r)=n;
    U = zeros(1,N);
    U = NumInhom(h,0,rmax,0,h,n,U);
    % determine alpha
    Z = 1 - (rmax^2+2*rmax+1)*exp(-rmax);
    alpha = (Z - U(N))/(N-1);
    for i=1:N
        U(i) = alpha*(i-1) + U(i);
    end
end