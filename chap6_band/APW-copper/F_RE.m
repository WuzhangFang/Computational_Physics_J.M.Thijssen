function F = F_RE(first,last,L,E,h)
    % RHS of the differential radial equation (6.27)
    F = zeros(1,last);
    for i=first:last
        r = i*h;
        F(i) = L*(L+1)/r^2 - 2*(V(r)/r + E);
    end % endfor
end