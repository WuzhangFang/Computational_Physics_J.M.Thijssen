% 6.7.4 Free particle in a box
L = 5;
b = 2 * pi / L;

% basis sets
K = [0 0 0; 
     b 0 0; 
     -b 0 0; 
     0 b 0; 
     0 -b 0; 
     0 0 b; 
     0 0 -b];
N = size(K,1);
H = zeros(N,N);
for i=1:7
    H(i,i) = norm(K(i,:))^2/2;
end

% diagonalize the Hamiltonian
[V, E] = eig(H);

disp(sum(E));

% calcualte the density (6.70)

Omega = L^3; % volume

% occupancy: 2 electrons at first energy level and 2 electrons at other
% levels
f = [2 1/3 1/3 1/3 1/3 1/3 1/3];

r = [1; 2; 1];

PhiR = getPhi(V, K, r);
PhiR2 = zeros(size(PhiR));

for i=1:size(PhiR,1)
    PhiR2(i) = 1/Omega * PhiR(i)' * PhiR(i);
end

n = f * PhiR2;

fprintf("The density is %f.\n", n);