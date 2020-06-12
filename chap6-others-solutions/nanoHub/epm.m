%Reference : Aaron Danner tutorial & Mathematica code (NUS/UIUC), & Kevin D. Welsher FORTRAN
%code
%Author: Muhanad Zaki
%Ref.(Danner): http://www.ece.nus.edu.sg/stfpage/eleadj/pseudopotential.htm
%Ref.(Welsher): http://large.stanford.edu/courses/ap272/welsher1/
%Empirical Pseudopotential Form factors are from Marvin L. Cohen and T.K. Bergstresser,
%“Band Structures and Pseudopotential Form Factors forFourteen Semiconductors of the Diamond and Zinc-blende Structures,”
% Phys. Rev. 141, 2, p. 556, 1966.

clear;
T=[1/8 1/8 1/8];   %Structure factor
%Start generating Brillouin zone vectors
n=(5^3)-1 ;
veck=zeros(3,125);
i=1;
for m=0:1:10  k(:,i)=[0.5-m/20 0.5-m/20 0.5-m/20] ;i=i+1; end;
for m=11:1:20 k(:,i)=[(m-10)/10 0 0] ; i=i+1;end;
for m=21:1:26 k(:,i)=[1 (m-20)/10 0] ; i=i+1;end;
for m=27:1:30 k(:,i)=[1-(m-25)/20 0.5+(m-25)/20 0] ; i=i+1;end;
for m=31:1:40 k(:,i)=[(0.75-(m-30)/(40/3)) (0.75-(m-30)/(40/3)) 0] ; i=i+1;end;

%List of constants
m0=9.11E-31; %free electron mass in Kg
a=5.43E-10;   %Temperature dependent
q=1.6E-19;
hbarJs=1.054E-34;
hbareV=6.581E-16;

% Start calculating the hamiltonian matrix elements,eigenvalues of this matrix 
% are the eigenenergies, E(k)

for d=1:1:41
    veck=k(:,d);
    for i=-62:1:62
        for j=-62:1:62
            if i==j ;
                m=i;
                G=Ggen(m);
                K(:,j+63)=G;
                DotProduct=dot((G(:)+veck(:)),(G(:)+veck(:))) ;
                A(i+63,j+63)=((hbarJs*hbareV)/(2*m0))*(abs(DotProduct*(2*pi/a)^2));

            else

                m=i-j;
                G=Ggen(m);
                K(:,j+63)=G;
                V=Vpseudo(G);
                A(i+63,j+63)=V;

            end;

        end;

    end;

    eval(['energy' num2str(d) '=eig(A)';]);

end;

%Save the programme output
for p=1:1:41

    eval(['e' num2str(p) '=energy' num2str(p) '(1:15)';]);
    eval(['Energies(:,p)' '=e' num2str(p) '(1:15)';]);
end;

% save 'C:\Energy.mat' Energies -ASCII      

for p=0:1:40
    kv=zeros(15);
    eval(['kv' num2str(p) '(1:15)' '=p';]);

end;

kays=zeros(15,41);
for n=0:1:40

    kays(:,n+1)=n;

end;
%Plot E(k) vs. k(0-40)
for p=1:1:41
    hold on;
    scatter(kays(:,p),Energies(:,p));
    grid on;
    title('bandstructure of Silicon');
    xlabel('k-vector');
    ylabel('Energy (eV)');
end;
hold off;


