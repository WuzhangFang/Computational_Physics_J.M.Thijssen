% use four Gaussian basis functions: e^(-alpha*r^2)
% where alpha has been optimized
alpha = [0.298073, 1.242567, 5.782948, 38.474970];

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
        A(i,j)=(-4*pi)/(alpha(i)+alpha(j));
    end
end

H = T + A;

% Hartree part
Q = zeros(N,N,N,N); % 4-dimensional array
for p=1:N
    for q=1:N
        for r=1:N
            for s=1:N
                sum_alpha = alpha(p)+alpha(q)+alpha(r)+alpha(s);
                Q(p,q,r,s) = (2*pi^2.5)/((alpha(p)+alpha(q))*(alpha(r)+alpha(s))*sqrt(sum_alpha)); % (4.17)
            end
        end
    end
end

% initialization
C=[1; 1; 1; 1]; % initial wavefunctions
C=normalize(C,S); % (4.19) normalize C
deltaE=1;
E=1;

% self-consistent calculations
% updating vector C until the change of energy is smaller than thershold.
while deltaE > 0.0001
   F=H+f(Q,C,N);
   [V,D] = eig(F,S); % (4.14)
   deltaE=abs(E-min(min(D)));
   E = min(min(D));
   [~,index]=find(D==E); % find the column corresponds to minimal E
   C=V(:,index);
   C=normalize(C,S);
end

% find ground state energy
E_G = ground(H,Q,C,N);

fprintf("The ground state energy of Hellium atom is %.8f Hartree.\n", E_G);

