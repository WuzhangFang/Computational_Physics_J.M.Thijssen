function C=normalize(C,S)
%normalize eigenvecors C based on the overlap matrix S
const = C'*S*C;
C = C/sqrt(const);
end