function y = F0(t)
    if t==0
        y=1;
    else
        y = t^(-1/2)*sqrt(pi)*erf(sqrt(t))/2;
    end
end