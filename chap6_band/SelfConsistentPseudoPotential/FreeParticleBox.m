% 6.7.4 Free particle in a box
L = 5;
b = 2 * pi / L;

% basis sets
K = [0 0 0; b 0 0; -b 0 0; 0 b 0; 0 -b 0; 0 0 b; 0 0 -b];

H = zeros(7,7);
for i=1:7
    H(i,i) = norm(K(i,:))^2/2;
end

E = eig(H);
