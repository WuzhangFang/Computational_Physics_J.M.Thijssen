function [G]=Ggen(m)
%This function generates one G vector
%G is the reciprocal lattice vector
%One G is generated per hamiltonian matrix element
%The function argument is only one number, m
%m=i if i==j, else m=i-j


lim=5;
n=lim^3-1;

h=floor((m+n/2)/(lim^2))-floor(lim/2);
k=floor(mod((m+n/2),lim^2)/lim)-floor(lim/2);
l=mod((m+n/2),lim)-floor(lim/2);

b1=[-1 1 1];  %basis vectors
b2=[1 -1 1];
b3=[1 1 -1];

B1(1)=h*b1(1);
B2(1)=k*b2(1);
B3(1)=l*b3(1);
B1(2)=h*b1(2);
B2(2)=k*b2(2);
B3(2)=l*b3(2);
B1(3)=h*b1(3);
B2(3)=k*b2(3);
B3(3)=l*b3(3);

G(1)=(B1(1)+B2(1)+B3(1)) ;
G(2)=(B1(2)+B2(2)+B3(2)) ;
G(3)=(B1(3)+B2(3)+B3(3)) ;

G=[G(1) G(2) G(3)];



