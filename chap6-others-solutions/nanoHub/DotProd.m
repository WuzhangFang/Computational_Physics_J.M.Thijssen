function [DotProd]=DotProd(G,veck)
%This function performs the dot product of
%the reciprocal lattice vector plus the brillouin zone vector


DotProd=dot(G+veck,G+veck);
