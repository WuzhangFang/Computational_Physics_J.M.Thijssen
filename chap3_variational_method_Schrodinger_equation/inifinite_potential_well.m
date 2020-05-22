
N = 5; % number of basis sets

% overlap matrix
S = zeros(N,N);
% initialization (3.17)
for m=0:N-1
    for n=0:N-1
        if rem(m+n,2)==0
            S(m+1,n+1)=(2/(n+m+5))-(4/(n+m+3))+(2/(n+m+1));
        end
    end
end

% Hamiltonian matrix
H = zeros(N,N);
% initialization (3.18)
for m=0:N-1
    for n=0:N-1
        if rem(m+n,2)==0
            H(m+1,n+1)=(-8)*((1-m-n-2*m*n)/((m+n+3)*(m+n+1)*(m+n-1)));
        end
    end
end

% solving the generalized eigenvalue equation
[V,D] = eig(H,S);
display(D);

