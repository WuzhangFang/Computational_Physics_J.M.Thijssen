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
CalcPoints (X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5)
double X1, Y1, X2, Y2, *X3, *Y3, *X4, *Y4, *X5, *Y5;
{
   *X3 = (2*X1 + X2)/3;
   *Y3 = (2*Y1 + Y2)/3;
   *X4 = (X1 + 2*X2)/3;
   *Y4 = (Y1 + 2*Y2)/3;
   *X5 = 0.5*(X1+X2-sqrt(3.0)/3*(Y2-Y1));
   *Y5 = 0.5*(Y1+Y2-sqrt(3.0)/3*(X1-X2));
}


void 
Koch(X1, Y1, X2, Y2)
double X1, Y1, X2, Y2;
{
  double Length, X3, Y3, X4, Y4, X5, Y5; 

  CalcPoints(X1, Y1, X2, Y2, &X3, &Y3, &X4, &Y4, &X5, &Y5);
  Length = sqrt((X2-X1)*(X2-X1)+ (Y2-Y1)*(Y2-Y1));
  
  if (Length>=0.05) 
    {
       Koch(X1, Y1, X3, Y3);
       Koch(X3, Y3, X5, Y5);
       Koch(X5, Y5, X4, Y4);
       Koch(X4, Y4, X2, Y2);
    }
  else
    Draw(X1, Y1, X2, Y2);
}


void 
main()
{
  PI = 4.0*atan (1.0);
  HalfSqrt3 = 0.5*sqrt(3.0);
  Sqrt2 = sqrt(2.0);
  Phi = atan(4.0+2.0*HalfSqrt3)-0.25*PI;
  TopScale = sqrt(1.25+HalfSqrt3);
  InitPlot ("lightblue", 800,400, "koch.ps", 2);
  PutStartButton();
  PutStopButton();
  Framing (-7.0, -1.0, 7.0, 6.0);
  Koch(-5.0, 0.0, 5.0, 0.0);
  SetFont("Times", "medium", "i", 18);
  MoveTo (3.0, 1.0);
  WriteText("Dit is de Koch-kromme");
  EndPlot();
}
