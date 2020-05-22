function U = NumInhom(Delta,StartR,MaxR,PhiStart,PhiNext,RHS,U)
    % Solve U''(r)=RHS(r) 
    % Delta: integration step
    % StartR and MaxR: integration boundary
    % PhiStart and PhiNext: first two intial points
    MaxI = (MaxR-StartR)/Delta;
    DeltaSq=Delta*Delta;
    Fac = DeltaSq/12;
    WPrev=PhiStart-Fac*RHS(1);
    U(1)=PhiStart;
    Phi=PhiNext;
    U(2)=PhiNext;
    W=Phi-Fac*RHS(2);
    for I=2:MaxI-1
         WNext = W*2 - WPrev + DeltaSq*RHS(I);
         WPrev = W;
         W     = WNext;
         Phi   = W+Fac*RHS(I+1);
         U(I+1) = Phi;
    end
    
end