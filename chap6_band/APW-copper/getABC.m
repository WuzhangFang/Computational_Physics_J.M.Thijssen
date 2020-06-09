function [A,B,C] = getABC(k,a,Lmax,R,K)
    % get the A, B, C matrix (6.35)
    % K: NKx3
    Vol = 3*a^3/4;
    factor = 2*pi*R^2/Vol;
    NK = size(K,1); % get the number of K vectors
    A = zeros(NK,NK);
    B = zeros(NK,NK);
    C = zeros(NK,NK,Lmax+1);
    for i=1:NK
        for j=1:NK
            Ki=K(i,:);
            Kj=K(j,:);
            qi=Ki+k;
            qj=Kj+k;
            di=norm(qi);
            dj=norm(qj);
            dij=norm(Ki-Kj);
            if dij<1E-5
                A(i,j)=-2*factor*R/3 + delta(i,j);
            else
                A(i,j) = -2*factor*Bessel(1,dij*R)/dij + delta(i,j);
            end
            B(i,j) = 0.5*A(i,j)*dot(qi,qj);
            for L=1:Lmax+1
                C(i,j,L)=(2*(L-1)+1)*factor*legendreP(L-1,dot(qi,qj)/(di*dj)) ...
                        *Bessel(L-1,di*R)*Bessel(L-1,dj*R);
            end %L    
         end %j
    end %i

end