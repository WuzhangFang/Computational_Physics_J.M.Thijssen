function solution = Verlet(h, iStart, iEnd, F, uStart, uNext, solution)
% u''(r)=F(r)
% h: integration step
% start: index of the starting point
% end: index of the ending point
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

uPrev = uStart;
solution(iStart) = uStart;

u = uNext;
solution(iStart+iStep) = uNext;

for i=iStart+iStep:iStep:iEnd-iStep
    % iStep is the central index
    uNext = 2*u - uPrev + h2*F(i)*u; % (2.12)
    uPrev = u;
    u = uNext;
    solution(i+iStep) = u;
end


end