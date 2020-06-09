function E=eigk(k,a,Lmax,R,K)
    % caculate the energy bands at k
    NK=size(K,1);
    [A,B,C] = getABC(k,a,Lmax,R,K);
    pDet = 0;
    i = 0;
    steps = 1;
    estep = 0.38/steps;
    for j=1:steps
        energy = -0.2 + (j-1)*estep;
        Hmat = H(energy,A,B,C,R,NK);
        Det = det(Hmat-energy*eyes(NK));
        if Det*pDet < 0
            i=i+1;
            E(i)=energy-Det*estep/(Det-pDet);
        end
        pDet=Det;
    end %for  
       
end