function [u,uprime] = atom(E, L, R)
    % solve the radial eqution (6.27)
    % return u(R) and u'(R)
    N = 1000/(L+1);
    h = R/N;
    % get array of RHS
    F = F_RE(1,N,L,E,h);
    % initial values
    if L==0
        uStart = 2*h^2*29/12;
    else
        uStart = 0;
    end %if
    
    uNext = h^(L+1);
    % call Numerov
    s = zeros(1,N);
    s = Numerov(h,0,N,F,uStart,uNext,s);
    u = s(N)/R;
    uprime = (s(N)-s(N-1))/h + 0.125*h*(3*F(N)*s(N)+F(N-1)*s(N-1));
    uprime = (uprime-u)/R;  
   
end