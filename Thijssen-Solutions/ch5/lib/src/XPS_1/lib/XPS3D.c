#include <math.h>
#include <stdio.h>

typedef double *dblptr;

static double Normal[3], v1[3], v2[3];
/* Normal is the "birds-eye" line.
   v1 and v2 are the vectors on which to project a 3d vector in order 
   to find the 2d x- and y-component
*/

double InnerProduct(a, b)
double a[3], b[3];
{
  double IP;
  IP = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return(IP);
}

void
CalcProject(p, d)
double p[3], d[3];
/* p is an arbitrary vector. It is projected onto the plane perpendicular to
  the bird's-eye axis */
{
  double Lambda;
  short int i;

  Lambda = InnerProduct(p, Normal);
  for (i=0; i<3; i++)
    d[i] =p[i]-Lambda*Normal[i];
}

double
NormOf(a)
double a[3];
{
  double Norm;
  Norm = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
  Norm = sqrt(Norm);
  return(Norm);
}

void
Normalise (a)
double a[3];
{
  double Norm;
  short int i;
  Norm = NormOf(a);
  Norm = 1.0/Norm;
  for (i=0; i<3; i++)
    a[i] = a[i]*Norm;
}

void
CalcVVecs()
{
  double zvec[3];
  zvec[0] = 0.0;
  zvec[1] = 0.0;
  zvec[2] = 1.0;
  CalcProject(zvec, v2);
  Normalise(v2);
  v1[2] = 0.0;
  v1[0] = -(v2[1]*Normal[2]-v2[2]*Normal[1]);
  v1[1] = (v2[0]*Normal[2]-v2[2]*Normal[0]);
  Normalise(v1);
  printf("%lf %lf %lf \n", v1[0], v1[1], v1[2]);
  printf("%lf %lf %lf \n", v2[0], v2[1], v2[2]);
  printf("%lf %lf %lf \n", Normal[0], Normal[1], Normal[2]);
}

void
Projection(D3Vec, D2Vec)
double D3Vec[3], D2Vec[2];
{
  double ProjectVec[3];
  CalcProject(D3Vec, ProjectVec);
  D2Vec[0] = InnerProduct(ProjectVec, v1);
  D2Vec[1] = InnerProduct(ProjectVec, v2);
}


void
InitBirdsEye(n1, n2, n3)
double n1, n2, n3;
{
  Normal[0] = n1;
  Normal[1] = n2;
  Normal[2] = n3;
  Normalise (Normal);
  CalcVVecs();
}

void
D3DrawLine (Sx, Sy, Sz, Ex, Ey, Ez)
double Sx, Sy, Sz, Ex, Ey, Ez;
{
  double Vec1[3], Vec2[3], D2V1[2], D2V2[2];
  Vec1[0] = Sx;
  Vec1[1] = Sy;
  Vec1[2] = Sz;
  Vec2[0] = Ex;
  Vec2[1] = Ey;
  Vec2[2] = Ez;
  Projection (Vec1, D2V1);
  Projection (Vec2, D2V2);
  Draw (D2V1[0], D2V1[1], D2V2[0], D2V2[1]);
}

/*
void
FillPoly3(points, npoints)
double points[];
int npoints;
{
   double *NewPoints, Vec[3], D2Vec[2];
   int i, j;
   NewPoints = (double *) calloc (npoints*2, sizeof(double));
   for (i=0; i<npoints; i++)
   {
     Vec[0] = points(3*i);
     Vec[1] = points(3*i+1);
     Vec[2] = points(3*i+2);
     Projection(Vec, D2Vec);
     NewPoints[2*i] = Vec[0];
     NewPoints[2*i+1] = Vec[1];
   }
   FillPolygon(NewPoints, npoints);
   cfree(NewPoints);
}
*/
   


void
FillVect(X, Y, Z, Vec3D)
double X, Y, Z, Vec3D[3];
{
  Vec3D[0] = X;
  Vec3D[1] = Y;
  Vec3D[2] = Z;
}


void
D3SetPoint (X, Y, Z)
double X, Y, Z;
{
  double Vec[3], D2Vec[2];
  FillVect(X, Y, Z, Vec);
  Projection(Vec, D2Vec);
  SetPoint(D2Vec[0], D2Vec[1]);
}
void
D3MoveTo (X, Y, Z)
double X, Y, Z;
{
  double Vec[3], D2Vec[2];
  FillVect(X, Y, Z, Vec);
  Projection(Vec, D2Vec);
  MoveTo(D2Vec[0], D2Vec[1]);
}

void
D3DrawTo (X, Y, Z)
double X, Y, Z;
{
  double Vec[3], D2Vec[2];
  FillVect(X, Y, Z, Vec);
  Projection(Vec, D2Vec);
  DrawTo(D2Vec[0], D2Vec[1]);
}


void
D3FillRect(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4)
double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4;
{
  double Vec1[3], Vec2[3], Vec3[3], Vec4[3], 
         D2Vec[2], D2Points[8];
  FillVect(X1, Y1, Z1, Vec1);
  FillVect(X2, Y2, Z2, Vec2);
  FillVect(X3, Y3, Z3, Vec3);
  FillVect(X4, Y4, Z4, Vec4);
  Projection(Vec1, D2Vec);
  D2Points[0] = D2Vec[0];
  D2Points[1] = D2Vec[1];
  Projection(Vec2, D2Vec);
  D2Points[2] = D2Vec[0];
  D2Points[3] = D2Vec[1];
  Projection(Vec3, D2Vec);
  D2Points[4] = D2Vec[0];
  D2Points[5] = D2Vec[1];
  Projection(Vec4, D2Vec);
  D2Points[6] = D2Vec[0];
  D2Points[7] = D2Vec[1];
  FillPolygon(D2Points, 4);
}

/*

void
DrawF(Xll, Yll, Xur, Yur, Nx, Ny)
double Xll, Yll, Xur, Yur;
int Nx, Ny;
{
  double XB, XF, YB, YF, hx, hy, x1, x2, y1, y2, z1, z2, z3, z4;
  int i, j;
  if ((Xll>Xur) || (Yll>Yur))
  {
    printf("Sorry, Xll must be < Xur and Yll < Yur \n");
    exit(0);
  }
  if (Normal[2]>0)
  {
    XB = Xll;
    YB = Yll;
    XF = Xur;
    YF = Yur;
  }
  hx = (XF-XB)/Nx;
  hy = (YF-YB)/Ny; 
  printf("%lf %lf\n", hx, hy);
  for (i=0; i<Nx; i++)
  {
    x1 = XB+i*(XF-XB)/Nx;
    x2 = x1+hx;
    for (j=0; j<Ny; j++)
    {
      y1= YB+j*(YF-YB)/Ny;
      y2= y1+hy;
      z1 = f(x1, y1);
      z2 = f(x2, y1);
      z3 = f(x1, y2);
      z4 = f(x2, y2);
      SetFastColor(3);
      D3FillRect(x1, y1, z1, x2, y1, z2, 
                x2, y2, z4, x1, y2, z3); 
      SetFastColor(6);
      D3DrawLine (x2, y2, z4, x1, y2, z3);
      D3DrawLine (x2, y2, z4, x2, y1, z2);
      D3DrawLine (x1, y1, z1, x1, y2, z3);
      D3DrawLine (x1, y1, z1, x2, y1, z2);
    }
  }
}
*/
  
void
D3DrawSurf(Psi, Size, MaxSize)
dblptr *Psi;
int Size;
{
  double x1, y1, x2, y2, z1, z2, z3, z4, Fac;
  int i, j;
  
  Fac = 32.0/((double)Size);
  for (i=0; i<Size-1; i++)
  {
    x1 = Fac*i;
    x2 = x1 + Fac;
    for (j=0; j<Size-1; j++)
    {
      y1 = Fac*j; 
      y2 = y1 + Fac;
      z1 = Psi[i][j];
      z2 = Psi[i+1][j];
      z3 = Psi[i][j+1];
      z4 = Psi[i+1][j+1];
      SetFastColor(3);
      D3FillRect(x1, y1, z1, x2, y1, z2, 
                 x2, y2, z4, x1, y2, z3); 
      SetFastColor(6);
      D3DrawLine (x2, y2, z4, x1, y2, z3);
      D3DrawLine (x2, y2, z4, x2, y1, z2);
      D3DrawLine (x1, y1, z1, x1, y2, z3);
      D3DrawLine (x1, y1, z1, x2, y1, z2);
    }
  }
}

