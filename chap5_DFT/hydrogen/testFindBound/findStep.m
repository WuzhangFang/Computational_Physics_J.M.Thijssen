function low = findStep(low,step,f)
    fPrev = f(low);
    while (f(low+step)*fPrev) > 0
        low=low+step;
    end
end