function E = interpolate(low,high,prec,f,V,h,N)
    leftX = low;
    rightX = high;
    leftY = f(low,V,h,N);
    rightY = f(high,V,h,N);
    midY = 1;
    while abs(midY) > prec
        midX = (rightX*leftY-leftX*rightY)/(leftY-rightY); % (A.6)
        midY = f(midX,V,h,N);
        if midY*leftY>0
            leftY = midY;
            leftX = midX;
        else
            rightY = midY;
            rightX = midX;
        end % endif
    end % endwhile
    E = midX;
end

