function [E, X, dist] = band(k1,k2,Nk,a,Kvectors)
    % get the band between k1 and k2
    NK = size(Kvectors,1); % get the size of Kvectors
    E = zeros(Nk,NK);
    for nk=1:Nk
        k=(k2-k1)*(nk-1)/(Nk-1)+k1; % get the k vector
        disp(k);
        E(nk,:)=eigk(k,a,Kvectors);
    end
    dist=norm(k2-k1); % get the distance 
    X=linspace(0,dist,Nk);
end