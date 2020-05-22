function Etot = getTotalEnergy(Eband, U, n, N,h)
    % Eband: band energy 
    % U: Hartree potential
    % n: density
    % r: grids of radius
    Ehar=calInt(-U.*n,N,h); % (5.82)
    Etot = 2*Eband-Ehar;
    
end