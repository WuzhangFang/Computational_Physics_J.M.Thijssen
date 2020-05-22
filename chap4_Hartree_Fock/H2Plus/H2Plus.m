% Program to calculate the ground state energy of H2Plus
% use eight Gaussian basis functions: e^(-alpha*r^2)
% first four are around nucleus A and the last four are around nucleus B

alphaA = [13.00773, 1.962079, 0.444529, 0.1219492];
alphaB = [13.00773, 1.962079, 0.444529, 0.1219492];
alpha = [alphaA alphaB];

% number of basis fucntions
N = length(alpha);

% initialize the position of nucleus A and B
Ra = 2; % can be any value
d = 1; % distance between nucleus A and B
Rb = Ra + d;

R = [Ra, Ra, Ra, Ra, Rb, Rb, Rb, Rb]; % index corresponding to basis functions

% overlap matrix
S = zeros(N,N);
for i=1:N
    for j=1:N
        S(i,j)=S_pq(alpha(i),R(i),alpha(j),R(j));
    end
end

% kinetic matrix
T = zeros(N,N);
for i=1:N
    for j=1:N
        T(i,j)=T_pq(alpha(i),R(i),alpha(j),R(j));
    end
end

% Colomb matrix
% Should be symmetric
A = zeros(N,N);
for i=1:N
    for j=1:N
        A(i,j) = A_pq(alpha(i),R(i),alpha(j),R(j),Ra) + A_pq(alpha(i),R(i),alpha(j),R(j),Rb);
    end
end

H = T + A;

[V,D] = eig(H,S);
fprintf("The ground state energy is %.6f Hartree.\n", min(min(D)));

