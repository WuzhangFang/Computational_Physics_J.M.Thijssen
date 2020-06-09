% caculate the band structure of FCC copper using APW

% muffin tin radius
R = 2.41191;
% lattice constant
a = 6.822;
Lmax=3;
% cutoff for basis set
% for Kcutoff=3, size(Kvec,1)=27
% for Kcutoff=4.5, size(Kvec,1)=113
Kcutoff = 3; % cutoff for basis set, for Kcutoff=3, size(K)=27
b = 2*pi/a;
% high-symmetry points
G = b*[0 0 0];
X = b*[1 0 0];
K = b*[3/4 3/4 0];
W = b*[1 1/2 0];
L = b*[1/2 1/2 1/2];
% number of k-points
Nk = 1;
Kvec = getKvectors(a,Kcutoff); 
% caculate the band structure between two high-symmetry points
[E,X,d] = band(G,X,Nk,a,Lmax,R,Kvec);
% plotting
plot(X,E);



