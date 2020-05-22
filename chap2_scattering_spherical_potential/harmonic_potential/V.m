%harmonic potential
function v = V(r)
rmax = 1; 
if r <= rmax
    v = r ^2;
else
    v =0;
end
