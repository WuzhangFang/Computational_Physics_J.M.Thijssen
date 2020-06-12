#ifndef _KMESH_
#define _KMESH_
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <algorithm>
#include <complex>
#include <map>
#include <algorithm>

using namespace std;

// Three dimensional vector is derived from usual vector<double> class
// It defined some convenient functions for 3D vectors.
class dvector3 : public vector<double>
{
public:
  dvector3(double x, double y, double z) : vector<double>(3)
  { (*this)[0]=x; (*this)[1]=y;  (*this)[2]=z; }
  dvector3() : vector<double>(3){};
  double operator*(const dvector3& a) const {return (*this)[0]*a[0]+(*this)[1]*a[1]+(*this)[2]*a[2];}
  double length() const{return sqrt(sqr((*this)[0])+sqr((*this)[1])+sqr((*this)[2]));}
  bool operator==(const dvector3& a){return ((fabs((*this)[0]-a[0])<1e-10) && (fabs((*this)[1]-a[1])<1e-10) && (fabs((*this)[2]-a[2])<1e-10));}
  friend dvector3 operator - (const dvector3& a, const dvector3& b);
  friend dvector3 operator + (const dvector3& a, const dvector3& b);
};
inline dvector3 operator*(const dvector3& a, double x)
{ return dvector3(a[0]*x,a[1]*x,a[2]*x);}
inline dvector3 operator*(double x, const dvector3& a)
{ return dvector3(a[0]*x,a[1]*x,a[2]*x);}
inline dvector3 operator/(const dvector3& a, double x)
{ return dvector3(a[0]/x,a[1]/x,a[2]/x);}

inline std::ostream& operator<<(std::ostream& stream, const dvector3& a){
  int width = stream.width();
  stream << std::setw(width)<< a[0] << " " << std::setw(width) << a[1] << " " << std::setw(width) <<a[2]<<" ";
  return stream;
}
inline dvector3 cross(const dvector3& a, const dvector3& b)
{
  return dvector3(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]);
}
inline double Vproduct(const dvector3& a, const dvector3& b, const dvector3& c)
{
  return cross(a,b)*c;
}
inline dvector3 operator - (const dvector3& a, const dvector3& b)
{
  return dvector3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}
inline dvector3 operator + (const dvector3& a, const dvector3& b)
{
  return dvector3(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}


int cmp(const dvector3& a, const dvector3& b)
{ return a.length()<b.length();}


void ChooseIrreducible(list<dvector3>& kp, vector<dvector3>& irkp, vector<double>& wkp)
{ // The function performs all symmetry operations of a cubic point-group to each k-point and
  // keeps only thos k-points which can not be obtained from another k-point by group operation.
  // These k-points are obviously irreducible.
  int kpsize=kp.size();// Number of all k-points
  list<dvector3> irkp0;// temporary list where irreducible k points will be stored
  list<int> wkp0;      // temporary weights
  while(kp.size()>0){  // continues until all k-points are grouped into irreducible classes 
    dvector3 tk(*kp.begin()); // we concentrate on the k-point which is the first in the list
    irkp0.push_back(tk);      // the first can be stored as irreducible
    wkp0.push_back(0);        // and the weights for this irreducible k-point is set to zero
    int& w0 = wkp0.back();    // reference to the weight is set for easier changing of the weight
    // We go over 48 symmetry operations of cubic system:
    // Each wector component can change sign: 2^3=8 possibilities
    // All permutations of components: 3!=6
    // Since the operations are independent, we have 3!*2^3=48 operations == number of cubic point group operations
    for (int ix=1; ix>=-1; ix-=2){  // three loops for all possible sign changes 
      for (int iy=1; iy>=-1; iy-=2){
	for (int iz=1; iz>=-1; iz-=2){
	  dvector3 nk(ix*tk[0], iy*tk[1], iz*tk[2]);
	  sort(nk.begin(),nk.end());// sorted so that next_permutation goes over all permutations
	  do{// all permutations of a vector are tried out
	    list<dvector3>::iterator iv = find(kp.begin(),kp.end(),nk);
	    // This permutation and sign change leads to some element still in the list of k-points?
	    if (iv!=kp.end()){
	      w0++;         // If yes, increase weight of irreducible k-point
	      kp.erase(iv); // and remove it from the list of the reducible k-points
	    }
	  } while(next_permutation(nk.begin(),nk.end()));// try next permutation
	}
      }
    }
  }
  // irreducible k-points are stored in the output vectors
  irkp.resize(irkp0.size());
  wkp.resize(wkp0.size());
  int j=0;
  for (list<dvector3>::const_iterator ik=irkp0.begin(); ik!=irkp0.end(); ik++) irkp[j++]=*ik;
  j=0;
  for (list<int>::const_iterator iw=wkp0.begin(); iw!=wkp0.end(); iw++) wkp[j++]=(*iw)/static_cast<double>(kpsize);
}

inline double kv0(int iq, int q)
{return (iq-static_cast<int>((q+1.5)/2)+1)/static_cast<double>(q);}

void Generate_All_K_Points(int q, const dvector3& b0, const dvector3& b1, const dvector3& b2, list<dvector3>& kp, int type=0)
{
  if (type==0){
    for (int ip=0; ip<q; ip++){
      double p = kv0(ip,q);
      for (int ir=0; ir<q; ir++){
	double r = kv0(ir,q);
	for (int is=0; is<q; is++){
	  double s = kv0(is,q);
	  dvector3 k = b0*p+b1*r+b2*s;
	  kp.push_back(k);
	}
      }
    }
  }else{
    for (int ip=0; ip<q; ip++){
      double p = (2*ip-q+1.)/(2*q);
      for (int ir=0; ir<q; ir++){
	double r = (2*ir-q+1.)/(2*q);
	for (int is=0; is<q; is++){
	  double s = (2*is-q+1.)/(2*q);
	  dvector3 k = b0*p+b1*r+b2*s;
	  kp.push_back(k);
	}
      }
    }
  }
}

#endif //_KMESH_
