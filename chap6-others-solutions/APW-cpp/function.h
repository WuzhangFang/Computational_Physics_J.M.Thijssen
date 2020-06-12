#ifndef FUNCTION_
#define FUNCTION_
#include <iostream>
#include <algorithm>
#include "assert.h"
//////////////////////////// Hierarchy of classes: ///////////////////////////////////////
//                             
//                                     base clas  = function<>
//                   |                                                |
//               function1D<>                                     funProxy<>
//
//                                       function2D<funProxy>
//

template<class T> class function;
template<class T> class function1D;
template<class T> class funProxy;
template<class T> class function2D;

//********************************************************************************//
// Base class for two derived classes: function1D<T> and function2D<T>.	  	  //
// It is also used as a proxy class for function2D. Function2D<T> consists of	  //
// arrays of function<T> rather than functions1D<T>.				  //
// Memory is allocated in a fortran-like fashion for better performance.	  //
// Linear interpolation is implemented with the operator() that takes one	  //
// argument (class intpar).							  //
//********************************************************************************//
template<class T>
class function{
protected:
  T *f;
  int N0, N;
public:
  T& operator[](int i) {Assert(i<N,"Out of range in function[]"); return f[i];}
  const T& operator[](int i) const {Assert(i<N,"Out of range in function[]"); return f[i];}
  const T& last() const {return f[N-1];}
  int size() const { return N;}
  int fullsize() const { return N0;}
  function& operator+=(const function& m);
  function& operator*=(const T& m);
  T* MemPt(){return f;}
  const T* MemPt() const{return f;}
  function& operator=(const T& c);
protected:
  function() : f(NULL), N0(0), N(0) {};
  explicit function(int N_) : N0(N_), N(N_) {};
  ~function(){};
  function(const function&){};
  template<class U> friend class function2D;
  template <class U> friend U scalar_product(const function<U>& f1, const function<U>& f2);
};

//******************************************************************//
// One dimensional functions derived from function<T>. It has it's  //
// own constructors and destructors.				    //
//******************************************************************//
template <class T>
class function1D : public function<T>{
public:
  function1D(){};
  explicit function1D(int N_);
  ~function1D();
  function1D(const function1D& m);
  void resize(int N_);
  function1D& operator=(const function1D& m);
  function1D& operator=(const T& c) {function<T>::operator=(c); return *this;}
  void Product(const function2D<T>& A, const function<T>& x, double alpha=1., double beta=0.);
};

template <class T>
class funProxy : public function<T>{
public:
  void Initialize(int N_, T* f_);
  void ReInitialize(int N_, T* f_);
  void resize(int N_);
  funProxy& operator=(const function<T>& m);
  ~funProxy(){};
};

//**********************************************************************//
// Two dimentional function<T> derived from function<T>. It consists	//
// of an array of function<T> rather tham function1D<T>.		//
// Constructor calls operator new and aferwords placement new operator	//
// to allocate the whole memory in one single large peace. 		//
//**********************************************************************//
template<class T>
class function2D{
protected:  
  void *memory;
  T* data;
  funProxy<T> *f;
  int N0, Nd0, N, Nd;
public:
  function2D() : memory(NULL), N0(0), Nd0(0), N(0), Nd(0) {};
  function2D(int N_, int Nd_);
  ~function2D();
  funProxy<T>& operator[](int i) {Assert(i<N,"Out of range in function2D[]"); return f[i];}
  const funProxy<T>& operator[](int i) const {Assert(i<N,"Out of range in function2D[]"); return f[i];}
  const T& operator()(int i, int j) const {Assert(i<N && j<Nd,"Out of range in function2D(i,j)"); return f[i].f[j];}
  T& operator()(int i, int j) {Assert(i<N && j<Nd,"Out of range in function2D(i,j)"); return f[i].f[j];}
  T* MemPt() { return data;}
  const T* MemPt() const { return data;}
  const int size_N() const {return N;}
  const int size_Nd() const {return Nd;}
  const int fullsize_N() const {return N0;}
  const int fullsize_Nd() const {return Nd0;}
  const int lda() const {return Nd0;}
  void resize(int N_, int Nd_);
  
  function2D& operator=(const function2D& m);
  function2D& operator+=(double x);
  function2D& operator+=(const function2D& m);
  function2D& operator-=(double x);
  function2D& operator-=(const function2D& m);
  function2D& operator=(const T& u);
  function2D& operator*=(const T& x);
  
  void TensorProduct(const function<T>& x, const function<T>& y, double alpha);
  void Product(const function2D& A, const function2D& B, const T& alpha, const T& beta);
  void DotProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
  void TProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
  void MProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
  void SymmProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
  void TSymmProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
};

// function ////////////////////////////////////////////////////////////////
template<class T> 
inline function<T>& function<T>::operator+=(const function& m)
{
  T_LOG(if (N!=m.size()) cerr << "Functions not of equal length! Can't sum!" << std::endl;)
  for (int i=0; i<N; i++) f[i] += m[i];
  return *this;
}

template<class T>
inline function<T>& function<T>::operator*=(const T& m)
{
  for (int i=0; i<N; i++) f[i] *= m;
  return *this;
}

template <class T>
inline function<T>& function<T>::operator=(const T& c)
{
  T_LOG(if (N<=0) cerr << "Size of function is non positive! "<<N<<std::endl;)
  for (int i=0; i<N; i++) f[i] = c;
  return *this;
}

// function1D ////////////////////////////////////////////////////////////
template<class T>
inline function1D<T>::function1D(int N_) : function<T>(N_)
{
  f = new T[N_];
}

template<class T>
inline function1D<T>::~function1D()
{
  delete[] f;
  f = NULL;
}

template<class T>
inline void function1D<T>::resize(int n)
{
  if (n>N0){
    if (f) delete[] f;
    f = new T[n];
    N0=n;
  }
  N = n;
}

template<class T>
inline function1D<T>::function1D(const function1D& m)
{
  resize(m.N);
  std::copy(m.f,m.f+N,f);
}

template <class T>
inline function1D<T>& function1D<T>::operator=(const function1D<T>& m)
{
  resize(m.N);
  std::copy(m.f,m.f+N,f);
  return *this;
}

// funProxy ///////////////////////////////////////////////////////////////
template <class T>
inline void funProxy<T>::Initialize(int N_, T* f_)
{
  N = N0 = N_; f = f_;
}

template <class T>
inline void funProxy<T>::ReInitialize(int N_, T* f_)
{
  N = N_; f = f_;
}

template <class T>
inline void funProxy<T>::resize(int N_)
{
  if (N_>N0) std::cerr<<"Can't resize funProxy, to small funProxy!"<<std::endl;
  else N=N_;
}

template <class T>
inline funProxy<T>& funProxy<T>::operator=(const function<T>& m)
{
  resize(m.size());
  std::copy(m.MemPt(),m.MemPt()+N,f);
  return *this;
}

#define HPoffset 8
// function2D ////////////////////////////////////////////////////////////
template<class T>
function2D<T>::function2D(int N_, int Nd_) : N0(N_), Nd0(Nd_), N(N_), Nd(Nd_) 
{
  memory = operator new (sizeof(funProxy<T>)*N0+sizeof(T)*Nd0*N0+HPoffset);
  
  Assert(memory!=NULL,"Out of memory");
  
  f = new (memory) funProxy<T>[N0];
  
  int offset = sizeof(funProxy<T>)*N0+HPoffset;
  data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
  
  for (int i=0; i<N0; i++) f[i].Initialize(Nd0,data+i*Nd0);
}

template<class T>
function2D<T>::~function2D()
{
  for (int i=0; i<N0; i++){
    f[i].~funProxy();
  }
  operator delete(memory);
  memory = NULL;
}

template <class T>
inline function2D<T>& function2D<T>::operator=(const function2D& m)
{
  if (m.N<=N0 && m.Nd<=Nd0){
    N = m.N; Nd = m.Nd;
    for (int i=0; i<N; i++) memcpy(f[i].f, m.f[i].f, sizeof(T)*Nd);
  } else{
    int msize = sizeof(funProxy<T>)*m.N+sizeof(T)*m.Nd*m.N+HPoffset;
    operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    memcpy(memory, m.memory, msize);
    N = N0 = m.N; Nd = Nd0 = m.Nd;
    f = new (memory) funProxy<T>[N];
    int offset = sizeof(funProxy<T>)*N+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N; i++) f[i].Initialize(Nd, data+i*Nd);
  }
  return *this;
}

template <class T>
inline void function2D<T>::resize(int N_, int Nd_)
{
  if (N_>N0 || Nd_>Nd0){
    //    clog<<"Deleting function2D and resizing from "<<N0<<" "<<Nd0<<" to "<<N_<<" "<<Nd_<<std::endl;
    int msize = sizeof(funProxy<T>)*N_ +sizeof(T)*Nd_*N_+HPoffset;
    operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    N = N0 = N_; Nd = Nd0 = Nd_;
    f = new (memory) funProxy<T>[N];
    int offset = sizeof(funProxy<T>)*N+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N; i++) f[i].Initialize(Nd, data+i*Nd);
  } else{
    N = N_; Nd = Nd_;
  }
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(double x)
{
  if (N!=Nd || !Nd || !N) {
    std::cerr << "Can't add number to non-square matrix!" << std::endl;
    return *this;
  }
  for (int i=0; i<Nd; i++) f[i][i] += x;
  return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(const function2D& m)
{
  if (N!=m.N || Nd!=m.Nd) {
    std::cerr << "Can't sum different matrices!" << std::endl;
    return *this;
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] += m[i][j];
  
  return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(double x)
{
  if (N!=Nd || !N || !Nd) {
    std::cerr << "Can't add number to non-square matrix!" << std::endl;
    return *this;
  }
  for (int i=0; i<Nd; i++) f[i][i] -= x;
  return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(const function2D& m)
{
  if (N!=m.N || Nd!=m.Nd) {
    std::cerr << "Can't sum different matrices!" << std::endl;
    return *this;
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] -= m[i][j];
  
  return *this;
}
template <class T>
inline function2D<T>& function2D<T>::operator=(const T& u)
{
  for (int i=0; i<N; i++) for (int j=0; j<Nd; j++) f[i].f[j]=u;
  return *this;
}
template <class T>
inline function2D<T>& function2D<T>::operator*=(const T& x)
{
  for (int i=0; i<N; i++) for (int j=0; j<Nd; j++) f[i][j] *= x;
  return *this;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const function<T>& f)
{
  int width = stream.width(); 
  for (int i=0; i<f.size(); i++) stream<<i<<" "<<std::setw(width)<<f[i]<<std::endl;
  return stream;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const function2D<T>& f)
{
  int width = stream.width(); 
  for (int i=0; i<f.size_N(); i++){
    for (int j=0; j<f.size_Nd(); j++)
      stream<<std::setw(width)<<f[i][j]<<" ";
    stream<<std::endl;
  }
  return stream;
}

template <class T, class functor>
T accumulate(const function2D<T>& data, functor& f)
{
  T sum=0;
  for (int i=0; i<data.size_N(); i++)
    for (int j=0; j<data.size_Nd(); j++)
      sum += f(data(i,j));
  return sum;
}

// template <class T>
// inline T sqr(const T& x){return x*x;}
template <class T>
inline T identity(const T& x){return x;}

extern "C" void dgetrf_(int* n1, int* n2, double* a, int* lda, int* ipiv,int* info);

double Determinant(function2D<double>& A)
{
  if (A.size_Nd()!=A.size_N()) {std::cerr<<"Can't compute determinant of nonquadratic matrix!"<<std::endl; return 0;}
  int info;
  int n = A.size_N();
  int lda = A.fullsize_Nd();
  function1D<int> ipiv(n);
  dgetrf_(&n, &n, A.MemPt(), &lda, ipiv.MemPt(), &info);
  if (info) {std::cerr<<"LU factorization complains : "<<info<<std::endl; return 0;}
  double det = 1;
  for (int i=0; i<n; i++) det *= ((ipiv[i]==i) ? 1 : -1) * A(i,i);
  return det;
}

#endif
