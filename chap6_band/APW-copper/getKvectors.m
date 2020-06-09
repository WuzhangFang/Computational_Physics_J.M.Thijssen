function K = getKvectors(a,Kcutoff)
    b1=2*pi/a*[-1 1 1];
    b2=2*pi/a*[1 -1 1];
    b3=2*pi/a*[1 1 -1];
    K=[];
    kmax=6;
    for L=-kmax:kmax
        for M=-kmax:kmax
            for N=-kmax:kmax
                k = L*b1 + M*b2 + N*b3;
                if norm(k) < 2*pi/a*Kcutoff
                    K = [K;k];
                end %if
            end %for N
        end %for M
    end %for L
end %function