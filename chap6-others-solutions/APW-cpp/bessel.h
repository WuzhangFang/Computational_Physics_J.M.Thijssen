#ifndef _BESSEL_
#define _BESSEL_

double DoubleFactorial(int l)
{// Calculates l!!
  int start = (l%2==0) ? 2 : 1;
  double val=start;
  for (int i=start+2; i<=l; i+=2) val*=i;
  return val;
}

double bessel_j(int l, double x)
{// Gives spherical bessel function with upward recursion
  double j0 = fabs(x)>1e-8 ? sin(x)/x : 1;// The error is of the order of x^3=1e-24
  if (l<=0) return j0; // in this case, we do not need j1
  double j1 = fabs(x)>1e-8 ? j0/x-cos(x)/x : x/3.;// Asymptotic expression of small x
  if (fabs(x)<1e-20) return 0;
  double j2 = j1;
  for (int i=2; i<=l; i++){
    j2 = j1*(2*i-1)/x - j0;
    j0 = j1;
    j1 = j2;
  }
  return j2;
}
double Bessel_j(int l, double x)
{// Gives spherical bessel function with downward recursion
  if (x>l) return bessel_j(l, x); // For large x, upward recursion works and is faster
  if (x<1e-20) return 0;
  int lstart = l + static_cast<int>(sqrt(40.*l)/2.);// This is an estimate where we need to start the recursion
  double j2=0, j1=1;
  double j0, jl, x1=1/x;// 1/x is stored for performance reasons
  for (int i=lstart; i>=0; i--){
    j0 = (2*i+3.)*x1*j1 - j2;
    if (i==l) jl = j0;
    j2 = j1;
    j1 = j0;
  }
  double true_j0 = sin(x)/x;
  return jl * true_j0/j0; // renormalizing the results
}
double dbessel_j(int l, double x)
{
  if (l==0 && fabs(x)<1e-20) return 0;
  if (fabs(x)<1e-5) return l*pow(x,l-1)/DoubleFactorial(2*l+1);
  return l*bessel_j(l,x)/x-bessel_j(l+1,x);
}
void dlog_bessel_j(int l, double x, double& dlogj, double& jl)
{// calculates x*d/dx log(j_l(x)) and j_l(x)
  if (fabs(x)<1e-5) { dlogj = l; jl = pow(x,l)/DoubleFactorial(2*l+1); return; }
  jl = bessel_j(l,x);
  dlogj= l-x*bessel_j(l+1,x)/jl;
}

double Legendre(int l, double x)
{
  switch(l){
  case 0: return 1;
  case 1: return x;
  case 2: return 1.5*x*x-0.5;
  case 3: return x*(2.5*x*x-1.5);
  case 4: return 0.375*(1-10*x*x*(1-1.1666666666666667*x*x));
  default:
    double p0 = x*(2.5*x*x-1.5);
    double p1 = 0.375*(1-10*x*x*(1-1.1666666666666667*x*x));
    double p2=0;
    for (int i=5; i<=l; i++){
      p2 = ((2*i-1)*x*p1-(i-1)*p0)/i;
      p0=p1;
      p1=p2;
    }
    return p2;
  }
}
#endif//_BESSEL_
