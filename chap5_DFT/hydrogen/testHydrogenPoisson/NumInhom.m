function solution = NumInhom(h, startR, maxR, uStart, uNext, F, solution)
% u''(r)=F(r)
% h: integration step
% F: array of F
% uStart: initial value of u 
% uNext: next value of u
% solution: u

N = (maxR-startR)/h;
h2 = h^2;
factor = h2 / 12;

wPrev = uStart-factor*F(1);
solution(1) = uStart;

u = uNext;
solution(2) = uNext;
w = u-factor*F(2);

for i=2:N-1
    % i is the central index
    wNext = 2*w - wPrev + h2*u*F(i); % (2.12)
    wPrev = w;
    w = wNext;
    u = w+factor*F(i+1);
    solution(i+1) = u;
end


end