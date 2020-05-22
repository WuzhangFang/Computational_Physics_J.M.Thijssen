% Hartree Program to calculate the ground state energy of H2
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

% Hartree part
Q = zeros(N,N,N,N); % 4-dimensional array
for p=1:N
    for r=1:N
        for q=1:N
            for s=1:N
                a = alpha(p) + alpha(q);
                b = alpha(r) + alpha(s);
                ra = (alpha(p)*R(p)+alpha(q)*R(q)) / (alpha(p) + alpha(q));
                rb = (alpha(r)*R(r)+alpha(s)*R(s)) / (alpha(r) + alpha(s));
                t = (alpha(p) + alpha(q))*(alpha(r) + alpha(s))*(ra-rb)^2 / (alpha(p)+alpha(q)+alpha(r)+alpha(s));
                s_pq = S_pq(alpha(p),R(p),alpha(q),R(q));
                s_rs = S_pq(alpha(r),R(r),alpha(s),R(s));
                Q(p,r,q,s) = 2*sqrt(a*b/(pi*(a+b)))*s_pq*s_rs*F0(t);
            end
        end
    end
end

% initialization
C=[1; 1; 1; 1; 1; 1; 1; 1;]; % initial wavefunctions
% normalization (4.19)
C=normalize(C,S); 
deltaE=1;
E=1;

while deltaE > 1E-6
   F=H+f(Q,C,N);
   [V,D] = eig(F,S);
   deltaE=abs(E-min(min(D)));
   E = min(min(D));
   [~,index]=find(D==E); % find the column index of lowest E
   C=V(:,index);
   C=normalize(C,S);
end

% find ground state energy
E_G = ground(H,Q,C,N);

fprintf("The ground state energy of H2 molecule by Hartree method (including nuclear repulsion +1) is %.8f Hartree.\n", E_G + 1);

