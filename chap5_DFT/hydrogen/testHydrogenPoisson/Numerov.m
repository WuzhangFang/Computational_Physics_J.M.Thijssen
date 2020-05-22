function solution = Numerov(h, iStart, iEnd, F, uStart, uNext, solution)
% u''(r)=F(r)
% h: integration step
% iStart: index of the starting point
% iEnd: index of the ending point
% F: array of F
% uStart: initial value of u 
% uNext: next value of u
% solution: u

if h < 0
    iStep = -1;
else
    iStep = 1;
end

h2 = h^2;
factor = h2 / 12;

wPrev = (1-factor*F(iStart)) * uStart;
solution(iStart) = uStart;

u = uNext;
solution(iStart+iStep) = uNext;
w = (1-factor*F(iStart+iStep))*u;

for i=iStart+iStep:iStep:iEnd-iStep
    % i is the central index
    wNext = 2*w - wPrev + h2*u*F(i); % (2.12)
    wPrev = w;
    w = wNext;
    u = w / (1-factor*F(i+iStep));
    solution(i+iStep) = u;
end


end