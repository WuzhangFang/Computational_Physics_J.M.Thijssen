#include <math.h>
#include <stdio.h>

static double PI, HalfSqrt3, Sqrt2, Phi, TopScale;


void
DrawSquare(X1, Y1, X2, Y2, X3, Y3)
double X1, X2, Y1, Y2, *X3, *Y3;
{
   double X4, Y4;

   *X3 = 0.5*(X1+X2+Y1-Y2);
   X4 = 0.5*(X1+X2-Y1+Y2);
   *Y3 = 0.5*(Y1+Y2-X1+X2);
   Y4 = 0.5*(Y1+Y2+X1-X2);
   Draw(X1, Y1, *X3, *Y3);
   DrawTo (X2, Y2);
   DrawTo (X4, Y4);
   DrawTo (X1, Y1);
}


void
CalcPoints (X1, Y1, Length, Angle, X2, Y2, X3, Y3)
double X1, Y1, Length, Angle, *X2, *Y2, *X3, *Y3;
{
   Angle = Angle + PI*0.25;
   Length = Length*Sqrt2; 
   *X2 = X1 + Length*cos(Angle);
   *Y2 = Y1 + Length*sin(Angle);
   Angle = Angle+Phi;
   Length = Length/Sqrt2;
   *X3 = X1 + Length*TopScale*cos(Angle);
   *Y3 = Y1 + Length*TopScale*sin(Angle);
}


void 
Pyth(X1, Y1, Length, Angle)
double X1, Y1, Length, Angle;
{
  double X2, Y2, X3, Y3, X5, Y5; 

  CalcPoints(X1, Y1, Length, Angle, &X2, &Y2, &X3, &Y3);
  
  DrawSquare(X1, Y1, X2, Y2, &X5, &Y5);

  if (Length>=0.1) 
    {
       Pyth(X5, Y5, 0.5*Length, Angle+PI/3);
       Pyth(X3, Y3, HalfSqrt3*Length, Angle-PI/6);
    }
}


void 
main()
{
  PI = 4.0*atan (1.0);
  HalfSqrt3 = 0.5*sqrt(3.0);
  Sqrt2 = sqrt(2.0);
  Phi = atan(4.0+2.0*HalfSqrt3)-0.25*PI;
  TopScale = sqrt(1.25+HalfSqrt3);
  InitPlot ("lightblue", 400,400, "pyth.ps", 2);
  PutStartButton();
  PutStopButton();
  Framing (-7.0, -1.0, 7.0, 13.0);
  Pyth(-3.0, 0.0, 2.0, 0.0);
  EndPlot();
}
