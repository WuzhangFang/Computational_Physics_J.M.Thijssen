% test the secant method (A.7)

low = -3;
step = 0.1;
prec = 0.01;

low = findStep(low,step,@ur);
high = low + step;

Emin = interpolate(low,high,prec,@ur);

display(Emin);
