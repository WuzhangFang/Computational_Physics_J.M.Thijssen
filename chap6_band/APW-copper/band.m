function [E, X, dist] = band(k1,k2,Nk,a,Lmax,R,K)
    % get the band between k1 and k2
    for nk=1:Nk
        k=(k2-k1)*(nk-1)/(Nk-1)+k1; % get the k vector
        E(nk,:)=real(eigk(k,a,Lmax,R,K));
    end
    dist=norm(k2-k1); % get the distance 
    X=linspace(0,dist,Nk);
end