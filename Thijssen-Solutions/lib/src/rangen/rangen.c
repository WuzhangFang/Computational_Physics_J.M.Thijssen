/******************************************************************
C* This file contains some random generator functions based on    *
C* the shift register principle.                                  *
C*                                                                *
C* This program is written in C.                                  *
C*                                                                *
C* The function "InitRand(InitJ)" must be called before anything  *
C* else. This initialises the random generator with seed InitJ.   *
C* The function "RealRand()" yields a random double precision,    *
C* uniformly distributed between 0 and 1.                         *
C* The function "InitRand(Num)" yields a random integer with      *
C* uniform distribution between 0 and Num-1                       *
C* The function "ExpRand(R1, R2)" yields two double precision     *
C* reals according to a Gaussian distribution. This is based on   *
C* the Box-Mulller algorithm.                                     *
C******************************************************************/
#include <math.h>

#define maxint 2147483647
#define factor 16807
#undef PI
/* This makes the program platform-independent */

static double RealFac = 0.4656612e-9;
static double PI;


static int rn[256];
static int jrand, pcount;

void
ran()
{
   jrand=jrand*factor;  
   if (jrand<0) 
     { jrand = jrand+maxint+1; }
}         




void
initrand_(initj)
int *initj;
{
   int i;

   jrand = *initj;
   for (i=0; i<256; i++)
      {  
         ran();
         rn[i]  = jrand;
      }
   pcount = 0;
   PI = 4.0*atan(1.0);
}



void
initrand(initj)
int *initj;
{
   int i;

   jrand = *initj;
   for (i=0; i<256; i++)
      {  
         ran();
         rn[i]  = jrand;
      }
   pcount = 0;
   PI = 4.0*atan(1.0);
}




int 
randgen()
{
   int pos1, pos2, scount;

   pos1 = (pcount);
   pos2 = (pcount+147) % 251;
   rn[pcount] = rn[pos1]^rn[pos2];
   scount = pcount;
   pcount = (pcount+1) % 251;
/*   printf("%d, %d\n", scount, rn[scount]); */
   return rn[scount];
 }




int 
intrand_(num)
int num;

{
   int IntRand;
   
   IntRand = RealFac*num*randgen();
   return (IntRand);
}


int 
intrand(num)
int num;

{
   int IntRand;
   
   IntRand = RealFac*num*randgen();
   return (IntRand);
}


double
realrand_()

{ 
   double RealRand;

   RealRand = randgen()*RealFac;
   return (RealRand);
}


double
realrand()

{ 
   double RealRand;

   RealRand = randgen()*RealFac;
   return (RealRand);
}




void
exprand (r1, r2)
double *r1, *r2;
{
/***** Subroutine uses  "realrand" to generate random numbers ********
  **** according to a Gaussian distribution. Width of this  **********
  **** distribution is 1. r1 and r2 are double precision's ***********/

double X, Y, Norm, RealRand, Phi, realrand();

Norm = realrand();
Phi = realrand()*2.0*PI;
Norm = sqrt(-2.0*log(Norm));
*r1 = Norm*cos(Phi);
*r2 = Norm*sin(Phi);
}



void
exprand_ (r1, r2)
double *r1, *r2;
{
/***** Subroutine uses  "realrand" to generate random numbers ********
  **** according to a Gaussian distribution. Width of this  **********
  **** distribution is 1. r1 and r2 are double precision's ***********/

double X, Y, Norm, RealRand, Phi, realrand();

Norm = realrand();
Phi = realrand()*2.0*PI;
Norm = sqrt(-2.0*log(Norm));
*r1 = Norm*cos(Phi);
*r2 = Norm*sin(Phi);
}




