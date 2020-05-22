% solve the radial equation of hydrogen (5.72) 
% using Numerov algorithm
% need to find E, such that ur(E)=0 corresponding to bound state

f = @ur;
x=[-0.2 -1]; % the solution boundaries
Emin = fzero(f,x);

display(Emin);