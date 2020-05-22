function E = interpolate(low,high,prec,f)
    leftX = low;
    rightX = high;
    leftY = f(low);
    rightY = f(high);
    midY = 1;
    while abs(midY) > prec
        midX = (rightX*leftY-leftX*rightY)/(leftY-rightY); % (A.6)
        midY = f(midX);
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

