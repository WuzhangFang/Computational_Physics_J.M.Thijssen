function Emin = findBound(U,h,N)
    % For U=0, E=-1.999999596745113, agrees with Fortran code
    low = -3;
    step = 0.1;
    prec=1E-8;
    
    % search for change of sign
    low = findStep(low,step,@getWaveFunction,U,h,N);
    high = low + step;

    % find root by interpolation
    Emin = interpolate(low,high,prec,@getWaveFunction,U,h,N);

end
