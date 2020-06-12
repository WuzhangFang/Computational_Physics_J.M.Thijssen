#ifndef NUMEROV
#define NUMEROV

#include <vector>

template <class container>
void Numerov(const container& F, int Nmax, double dx, std::vector<double>& Solution)
{
  double h2 = dx*dx;
  double h12 = h2/12;
  
  double w0 = (1-h12*F[0])*Solution[0];
  double Fx = F[1];
  double w1 = (1-h12*Fx)*Solution[1];
  double Phi = Solution[1];
  
  double w2;
  for (int i=2; i<Nmax; i++){
    w2 = 2*w1 - w0 + h2*Phi*Fx;
    w0 = w1;
    w1 = w2;
    Fx = F[i];
    Phi = w2/(1-h12*Fx);
    Solution[i] = Phi;
  }
}
template <class container>
void NumerovInhom(const container& U, int Nmax, double dx, std::vector<double>& Solution)
{
  double h2 = dx*dx;
  double h12 = h2/12;
  
  double w0 = Solution[0]-h12*U[0];
  double Ux = U[1];
  double w1 = Solution[1]-h12*Ux;
  double Phi = Solution[1];
  
  double w2;
  for (int i=2; i<Nmax; i++){
    w2 = 2*w1 - w0 + h2*Ux;
    w0 = w1;
    w1 = w2;
    Ux = U[i];
    Phi = w2+h12*Ux;
    Solution[i] = Phi;
  }
}
template <class container>
void NumerovGen(const container& F, const container& U, int Nmax, double dx, std::vector<double>& Solution)
{
  double h2 = dx*dx;
  double h12 = h2/12;
  
  double w0 = Solution[0]*(1-h12*F[0])-h12*U[0];
  double Fx = F[1];
  double Ux = U[1];
  double w1 = Solution[1]*(1-h12*Fx)-h12*Ux;
  double Phi = Solution[1];
  
  double w2;
  for (int i=2; i<Nmax; i++){
    w2 = 2*w1 - w0 + h2*(Phi*Fx+Ux);
    w0 = w1;
    w1 = w2;
    Fx = F[i];
    Ux = U[i];
    Phi = (w2+h12*Ux)/(1-h12*Fx);
    Solution[i] = Phi;
  }
}
#endif//NUMEROV
