function [j, n] = Bessel(l, r)
% [j, n] = Bessel(l, r): return 1st and 2nd kind of spherical Bessel functions
j = sqrt(pi/(2*r))*besselj(l+0.5, r);
n = sqrt(pi/(2*r))*bessely(l+0.5, r);
end