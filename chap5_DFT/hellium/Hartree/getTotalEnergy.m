function Etot = getTotalEnergy(Eband, U, n, N,h)
    % Eband: band energy 
    % U: Hartree potential
    % n: density
    % r: grids of radius
    r=(0:N-1)*h;
    Ehar=trapz(r,-U.*n); % (5.82)
    Etot = 2*Eband-Ehar;
    
end