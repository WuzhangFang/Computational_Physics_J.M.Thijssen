#ifndef _INETGRATE
#define _INTEGRATE

template <class T, class container>
T integrate4(const container& F, double dh, int Np)
{ // Forth order integration routine
  static double coeff[4] = {17./48., 59./48.,43./48.,49./48.};
  T sum = coeff[0]*F[0]+coeff[1]*F[1]+coeff[2]*F[2]+coeff[3]*F[3];
  for (int i=4; i<Np-4; i++) sum += F[i];
  sum += coeff[0]*F[Np-1]+coeff[1]*F[Np-2]+coeff[2]*F[Np-3]+coeff[3]*F[Np-4];
  return sum*dh;
}

template <class T>
inline T sqr(const T& a){return a*a;}

#endif //_INTEGRATE
