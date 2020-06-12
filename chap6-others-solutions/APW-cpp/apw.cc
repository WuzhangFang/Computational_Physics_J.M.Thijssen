#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <list>
#include "numerov.h"
#include "integrate.h"
#include "kmesh.h"
#include "bessel.h"
#include "function.h"

extern "C" {
  void ludcmp_(double* a, int* n, int* np, int* indx, double* d);
};

using namespace std;

double Determinant1(function2D<double>& A)
{// C++ wrapper function for calculating Determinant
  int si = A.size_N();
  function1D<int> ipivot(si);
  double D;
  ludcmp_(A.MemPt(), &si, &si, ipivot.MemPt(), &D);
  for (int j=0; j<si; j++) D = D*A(j,j);
  return D;
}


class FccLattice{ // Class for storing reciprocal lattice
  double LatConst, Volume;
  dvector3 a0, a1, a2;    // Primitive vectors of fcc lattice
  dvector3 b0, b1, b2;    // Primitive vectors of reciprocal lattice
  dvector3 GammaPoint, LPoint, KPoint, XPoint, WPoint; // Special points in 1IRB
  vector<dvector3> Kmesh, kmesh;
public:
  double Vol(){return Volume;}
  int Ksize(){return Kmesh.size();}
  int ksize(){return kmesh.size();}
  const dvector3& K(int i){return Kmesh[i];} // can not be changed, only read
  const dvector3& k(int i){return kmesh[i];} // can not be changed, only read
  FccLattice(double LatConst_) : LatConst(LatConst_)
  {
    a0 = dvector3(0.5*LatConst,0.5*LatConst,0);
    a1 = dvector3(0.5*LatConst,0,0.5*LatConst);
    a2 = dvector3(0,0.5*LatConst,0.5*LatConst);
    Volume = fabs(Vproduct(a0,a1,a2));// Volume
    clog<<"Volume is "<<Volume<<endl;
    b0 = (2*M_PI/Volume)*cross(a1,a0);
    b1 = (2*M_PI/Volume)*cross(a0,a2);
    b2 = (2*M_PI/Volume)*cross(a2,a1);
    // Special points in Brillouin zone
    double brs = 2*M_PI/LatConst;
    GammaPoint = dvector3(0,0,0);
    LPoint = dvector3(0.5*brs,0.5*brs,0.5*brs);
    KPoint = dvector3(0.75*brs,0.75*brs,0);
    XPoint = dvector3(1*brs,0, 0);
    WPoint = dvector3(1*brs,0.5*brs,0);
  }
  void GenerateReciprocalVectors(int q, double CutOffK)
  {
    // Many reciprocal vectors are generated and later only the sortest are used
    list<dvector3> Kmesh0;
    for (int n=-q; n<q; n++){
      for (int l=-q; l<q; l++){
	for (int m=-q; m<q; m++){
	  Kmesh0.push_back(n*b0+l*b1+m*b2);
	}
      }
    }
    Kmesh0.sort(cmp); // Sorting according to the length of vector. Shortest will be kept
    int Ksize=0;
    for (list<dvector3>::const_iterator l=Kmesh0.begin(); l!=Kmesh0.end(); l++,Ksize++) if (l->length()>CutOffK) break;
    Kmesh.resize(Ksize);
    int j=0;
    for (list<dvector3>::const_iterator l=Kmesh0.begin(); l!=Kmesh0.end() && j<Ksize; l++,j++) Kmesh[j]=*l;
    clog<<"K-mesh size="<<Kmesh.size()<<endl;
  }
  void ChoosePointsInFBZ(int nkp){// Chooses the path in the 1BZ we will use
    kmesh.resize(nkp);
    int N0=kmesh.size()/4;
    for (int i=0; i<N0; i++) kmesh[i]      = GammaPoint + (XPoint-GammaPoint)*i/(N0-1.);
    for (int i=0; i<N0; i++) kmesh[N0+i]   = XPoint + (LPoint-XPoint)*i/(N0-1.);
    for (int i=0; i<N0; i++) kmesh[N0*2+i] = LPoint + (GammaPoint-LPoint)*i/(N0-1.);
    for (int i=0; i<N0; i++) kmesh[N0*3+i] = GammaPoint + (KPoint-GammaPoint)*i/(N0-1.);
  }
};

// This is the parametrization for the effective potential from ....
double VeffP(double R)
{
  return 29*exp(-2.3151241717834*pow(R,0.81266614122432)+ 2.1984250222603e-2*pow(R,4.2246376280056))
    -0.15595606773483*R-3.1350051440417e-3*R*R+5.1895222293006e-2*pow(R,3)-2.8027608685637e-2*pow(R,4);
}

class PartialWave{// Class for solving SCH equation
  int Z;
  vector<double> Rmesh;
  vector<double> rhs_MT, ur, temp; // For solving SCH equation
public:
  PartialWave(int N, double RMuffinTin, int Z_): Z(Z_), Rmesh(N), rhs_MT(N), ur(N), temp(N)
  {
    // Building equidistant radial mesh
    double dh = RMuffinTin/(N-1.);
    for (int i=0; i<N; i++) Rmesh[i] = i*dh;
    clog<<"RmuffinTin="<<Rmesh[Rmesh.size()-1]<<endl;
  }
  double startSol(int Z, int l, double r) // good choice for starting Numerov algorithm
  { return pow(r,l+1)*(1-Z*r/(l+1));}
  double LogDerivative(double Enu, int l)
  {
    rhs_MT[0]=0;
    for (int i=1; i<Rmesh.size(); i++){
      double Veff = -VeffP(Rmesh[i])/Rmesh[i] + 0.5*l*(l+1)/sqr(Rmesh[i]);
      rhs_MT[i] = 2*(Veff-Enu);
    }
    double dh = Rmesh[1]-Rmesh[0];
    ur[0]=0;
    ur[1]=startSol(Z, l, dh);
    // Solving radial SCH equation
    Numerov(rhs_MT, ur.size(), dh, ur);
    // Normalizing the result  - here it is not necessary because we only need log derivative
    for (int i=0; i<Rmesh.size(); i++) temp[i] = sqr(ur[i]);
    double norm = 1./sqrt(integrate4<double>(temp, dh, temp.size()));
    for (int i=0; i<ur.size(); i++) ur[i] *= norm;
    // Here we estimate the derivative at the Muffin-Tin boundary
    int N0 = ur.size()-1;
    double v1 = rhs_MT[N0]*ur[N0];
    double v0 = rhs_MT[N0-1]*ur[N0-1];
    double dudr  = (ur[N0]-ur[N0-1])/dh + 0.125*dh*(3*v1+v0);
    double RMuffinTin = Rmesh[N0];
    return RMuffinTin*dudr/ur[N0] -  1;
  }
};

int main(int argc, char *argv[], char *env[])
{
  int Z=29;                    // Number of electrons in the Cu atom
  double E_start = -0.2;      // Where to start searching for Energy
  double dE = 1e-3;            // Step in serching for bound states
  int nE = 500;                // Number of enrgy steps when searching for energy bands
  double LatConst = 6.8219117; // Lattic constant
  double RMuffinTin = 2.41191; // Muffin tin radius - Touching spheres
  int lMax=5;                  // Maximum l considered in calculation
  int N = 1001;                // Number of points in radial mesh
  int nkp = 40;                // Number of k-points in 1BZ
  double CutOffK=3.;           // Largest lengt of reciprocal vectors K (only shorter vec. are taken into account)
  
  int i=0;
  while (++i<argc){
    std::string str(argv[i]);
    if (str=="-dE" && i<argc-1) dE = atof(argv[++i]);
    if (str=="-h" || str=="--help"){
      std::clog<<"**************** APW program for Cu-fcc **************\n";
      std::clog<<"**                                                  **\n";
      std::clog<<"**      Copyright Kristjan Haule, 18.10.2005        **\n";
      std::clog<<"******************************************************\n";
      std::clog<<"\n";
      std::clog<<"apwCu [-dE double] [] []\n" ;
      std::clog<<"Options:   -Z          Number of electrons ("<<Z<<")\n";
      std::clog<<"           -dE         Step in searching for states ("<<dE<<")\n";
      std::clog<<"*****************************************************\n";
      return 0;
    }
  }

    
  clog.precision(10);
  // For solving SCH equation
  PartialWave wave(N, RMuffinTin, Z);
  // Generates and stores momentum points
  FccLattice fcc(LatConst);                  // Information about lattice
  fcc.GenerateReciprocalVectors(4, CutOffK); // Reciprocal bravais lattice is builded
  fcc.ChoosePointsInFBZ(nkp);                // Chooses the path in the 1BZ we will use
  
  // Storage for matrices
  function2D<double> Olap0(fcc.Ksize(),fcc.Ksize());
  function2D<double> Ham0(fcc.Ksize(),fcc.Ksize()), Hmat(fcc.Ksize(),fcc.Ksize());
  vector<function2D<double> > Ham1(lMax+1);
  for (int il=0; il<=lMax; il++) Ham1[il].resize(fcc.Ksize(),fcc.Ksize());

  vector<double> dlogPsi(lMax+1);
  function2D<double> jl0(fcc.Ksize(), lMax+1), Dl0(fcc.Ksize(),lMax+1);
  
  vector<list<double> > spagety(fcc.ksize());
  
  // Overlap in the interstitials can be calculated outside
  for (int i=0; i<fcc.Ksize(); i++){
    Olap0(i,i) = 1 - 4*M_PI*sqr(RMuffinTin)*RMuffinTin/(3.*fcc.Vol());
    for (int j=i+1; j<fcc.Ksize(); j++){
      double KKl = (fcc.K(i)-fcc.K(j)).length();
      Olap0(i,j) = -4*M_PI*sqr(RMuffinTin)*bessel_j(1,KKl*RMuffinTin)/(KKl*fcc.Vol());
      Olap0(j,i) = Olap0(i,j);
    }
  }
  
  // Main loop over k points
  for (int ik=0; ik<fcc.ksize(); ik++){
    dvector3 k = fcc.k(ik);
    clog<<"k="<<ik<<endl;
    
    // Bessel functions can be calculated only ones for each K-point
    for (int iK=0; iK<fcc.Ksize(); iK++)
      for (int il=0; il<=lMax; il++)
	dlog_bessel_j(il, (k+fcc.K(iK)).length()*RMuffinTin, Dl0(iK,il), jl0(iK,il)); 
    // Parts of the Hamiltonian matrix which do not depend on energy, are calculated
    for (int iK=0; iK<fcc.Ksize(); iK++){
      for (int jK=0; jK<=iK; jK++){
	dvector3 qi(k+fcc.K(iK));
	dvector3 qj(k+fcc.K(jK));
	
	Ham0(iK,jK) = (0.5*qi*qj)*Olap0(iK,jK);
	Ham0(jK,iK) = Ham0(iK,jK);
	
	double qi_len = qi.length();
	double qj_len = qj.length();
	double argv = (qi_len*qj_len==0) ? 1. : qi*qj/(qi_len*qj_len);
	
	double cc = (2*M_PI*RMuffinTin/fcc.Vol());
	for (int il=0; il<=lMax; il++){
	  Ham1[il](iK,jK) = cc*(2*il+1)*Legendre(il,argv)*jl0(iK,il)*jl0(jK,il);
	  Ham1[il](jK,iK) = Ham1[il](iK,jK);
	}
      }
    }
    // Here we scan energy interval and look for possible solutions
    double Det0 = 0, Det1; // New and old determinant
    for (int ie=0; ie<nE; ie++){
      
      double Enu = E_start + ie*dE; // current energy
      
      for (int l=0; l<=lMax; l++) dlogPsi[l]  = wave.LogDerivative(Enu,l);

      for (int iK=0; iK<fcc.Ksize(); iK++){
	for (int jK=0; jK<=iK; jK++){
	  double sum=0;
	  //	  clog<<"Dl="<<(0.5*Dl0(iK,0)+0.5*Dl0(jK,0))<<endl;
	  for (int il=0; il<=lMax; il++) sum += Ham1[il](iK,jK)*dlogPsi[il];//-Dl0(iK,il));//-0.5*Dl0(iK,il)-0.5*Dl0(jK,il));
	  Hmat(iK,jK) = Ham0(iK,jK) - Enu*Olap0(iK,jK) + sum;
	  Hmat(jK,iK) = Hmat(iK,jK);
	}
      }
      Det1 = Determinant1(Hmat);
      if (Det1*Det0<0){ // solution bracketed and linear interpolation used!
	double Ee = Enu-dE*Det1/(Det1-Det0); // linear interpolation
	spagety[ik].push_back(Ee);
	clog<<Ee<<" "<<endl;
      }
      Det0 = Det1;
    }
  }

  // Printing of results
  for (int i=0; i<spagety.size(); i++){
    cout<<setw(10)<<static_cast<double>(i)/(spagety.size()-1)<<" ";
    for (list<double>::iterator l = spagety[i].begin(); l!=spagety[i].end(); l++) cout<<setw(10)<<*l<<" ";
    cout<<endl;
  }
  return 0;

}
