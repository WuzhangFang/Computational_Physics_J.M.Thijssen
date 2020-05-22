% figure 2.4
% energy in meV, radius in rho=3.57 angstrom
clear;

h=0.01; % integration step
L=6;
rmax=5;
emesh=1000;

E=linspace(0.1,3.5,emesh);
sigma=zeros(size(E));

% sum over L
for l=0:L
    for i=1:emesh
        sigma(i) = sigma(i) + crossSection(l,E(i),h,rmax);
    end
end

%plotting
plot(E,sigma);
xlabel('E (meV)')
ylabel('Total cross section (\rho^2)')
axis([0,3.5,5,50])
box on
saveas(gcf,'Fig2.4.png')
