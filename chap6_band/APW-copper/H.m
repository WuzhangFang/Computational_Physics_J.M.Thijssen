function Hmat = H(E,Lmax,A,B,C,R,NK)
    Hmat = zeros(NK,NK);
    for i=1:NK
        for j=1:NK
            Hmat(i,j)=-E*A(i,j)+B(i,j);
            for L=1:Lmax+1
                [u,uprime]=atom(E,L-1,R);
                Hmat(i,j)=Hmat(i,j)+C(i,j,L)*uprime/u;
            end %L
         end %j     
    end %i
end
