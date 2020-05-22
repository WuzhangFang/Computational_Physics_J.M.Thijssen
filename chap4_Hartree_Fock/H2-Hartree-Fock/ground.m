function E_G = ground(H,Q,C,N)
% get the ground state energy based on the eigenvectors
E_G = 2 * C'*H*C;
F=0;

for p=1:N
    for q=1:N
        for r=1:N
            for s=1:N
                F=F+Q(p,q,r,s)*C(p)*C(q)*C(r)*C(s); % 4.21
            end
        end
    end
end

E_G = E_G + F;

end