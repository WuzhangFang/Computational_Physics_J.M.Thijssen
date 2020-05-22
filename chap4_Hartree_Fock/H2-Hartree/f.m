function F = f(Q,C,N)
    F = zeros(N,N);
    for p=1:N
        for r=1:N
            for q=1:N
                for s=1:N
                    F(p,q)=F(p,q)+Q(p,r,q,s)*C(r)*C(s);
                end
            end
        end
    end
            

end