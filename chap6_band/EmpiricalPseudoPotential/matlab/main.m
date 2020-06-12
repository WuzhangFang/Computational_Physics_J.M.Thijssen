% calculating band structure of Silicon using pseudopotential method
a = 5.43; % lattice constant in a.u.

Nk = 31; % number of kpoints between two high-symmetry point
% high-symmetry points in BZ of fcc
G = [0 0 0];
X = 0.5*[1 0 0];

load("Kvectors");
% get the band structure
[E,xx,dist] = band(G,X,Nk,a,Kvectors);

% plotting
plot(xx,E);
axis([0,dist,-0.1,0.6])