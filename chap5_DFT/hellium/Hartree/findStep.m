function low = findStep(low,step,f,V,h,N)
    fPrev = f(low,V,h,N);
    while (f(low+step,V,h,N)*fPrev) > 0
        low=low+step;
    end
end