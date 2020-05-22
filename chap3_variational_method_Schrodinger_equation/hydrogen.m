% use four Gaussian basis functions: e^(-alpha*r^2)

alpha = [13.00773, 1.962079, 0.444529, 0.1219492];

% number of basis fucntions
N = length(alpha);

% overlap matrix
S = zeros(N,N);
for i=1:N
    for j=1:N
        S(i,j)=(pi/(alpha(i)+alpha(j)))^1.5;
    end
end

% kinetic matrix
T = zeros(N,N);
for i=1:N
    for j=1:N
        T(i,j)=3*(alpha(i)*alpha(j))*pi^1.5/(alpha(i)+alpha(j))^2.5;
    end
end

% Colomb matrix
A = zeros(N,N);
for i=1:N
    for j=1:N
        A(i,j)=(-2*pi)/(alpha(i)+alpha(j));
    end
end

H = T + A;


[V,D] = eig(H,S);
fprintf("The ground state energy is %.5f Hartree.\n", min(min(D)));

