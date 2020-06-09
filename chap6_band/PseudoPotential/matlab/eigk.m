function E = eigk(k,a,Kvectors)
    % caculate the energy bands at k
    % generate K vecotors: basis sets (l,m,n)
    NK = size(Kvectors,1);
    H = zeros(NK,NK);
    b = 2*pi/a;
    % fcc reciprocal lattice
    b1 = 0.5*[-1 1 1];
    b2 = 0.5*[1 -1 1];
    b3 = 0.5*[1 1 -1];
    % fill in the H matrix (6.58)
    for i=1:NK
        for j=1:NK
            deltaK = Kvectors(i,:)-Kvectors(j,:);
            normDeltaK = norm(deltaK);
            if normDeltaK == sqrt(3)
                H(i,j) = H(i,j) - 0.1121 * cosSum(deltaK);
            elseif normDeltaK == sqrt(8)
                H(i,j) = H(i,j) + 0.0276 * cosSum(deltaK);
            elseif normDeltaK == sqrt(11)
                H(i,j) = H(i,j) + 0.0362 * cosSum(deltaK);
            end %if
        end %j
        % G = l*b1 + m * b2 + n * b3
        G = Kvectors(i,1)*b1+Kvectors(i,2)*b2+Kvectors(i,3)*b3;
        ktotal = k + G;
        H(i,i) = 1/2*b^2*(norm(ktotal))^2;
    end %i
    E = eig(H);
    
end