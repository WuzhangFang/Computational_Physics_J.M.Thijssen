
function [V]=Vpseudo(G)
%This function calculates the pseudopotential
%corresponding to a certain G vector
%This is only in the case of i~=j (not equal)
%The function argument is the reciprocal lattice vector G
%The Form Factors are in this file.

DotProduct=dot(G(:),G(:));
if     DotProduct==3  Vs=(-0.21)*13.6059;  %13.6 is the Rydberg energy in eV's
elseif DotProduct==8  Vs=0.04*13.6059;
elseif DotProduct==11 Vs=0.08*13.6059;
else                  Vs=0;
end;



T=[1/8 1/8 1/8];

V=Vs*cos(2*pi*dot(G(:),T(:)));


