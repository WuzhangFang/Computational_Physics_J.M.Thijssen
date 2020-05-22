% solve the ground state l=0 for hellium using DFT
% with Hartree potential of one electron
% U''(r)=n(r) is not solved accurately

% paramaters and initialization
precision=1E-7;
rmax=10;
h=0.01;
N=rmax/h;

% self-consistent calculations
Eprev=-10;
E=-0.3;
U = zeros(1,N);
i=1;
while abs(Eprev-E) > precision
    Eprev = E;
    E = findBound(U,h,N);
    disp(i);
    disp(E);
    [u0,~,n] = getWaveFunction(E,U,h,N);
    U = getHartreePotential(n,h,N);
    i=i+1;
end

Etot = getTotalEnergy(E,U,n,N,h);
fprintf("The ground state energy of He by DFT is %.6f Hartree.\n", Etot);