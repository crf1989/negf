#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <assert.h>
using namespace std;
using namespace Eigen;

const double PI = 3.141592653589793;
const double BOLTZMANN_CONSTANT = 8.6173324e-5;

template<typename T>
T square (T a)
{
  return a*a;
}

double fermi (double e, double u, double T)
{
  assert (T >= 0);
  if (T == 0)
    if (e == u)
      return 0.5;
    else 
      return (e > u ? 0 : 1);
  else
    return 1/(exp((e-u)/(BOLTZMANN_CONSTANT*T))+1);
}

complex<double> sfg_solve (double e, double t)
{
  complex<double> lambda1 = 0.5*(e/t-sqrt((complex<double>)(square(e/t)-4)));
  complex<double> lambda2 = 0.5*(e/t+sqrt((complex<double>)(square(e/t)-4)));
  if (abs(lambda1) <= 1)
    return lambda1*t;
  else
    return lambda2*t;
}

int main (int argc, char* argv[])
{
  int L = 10;
  int N = 1000000;
  double h = 1;
  double dw = 4.1/N;
  double* T = new double[N];

  double TL = 500;
  double TR = 500;
  double uL = atof(argv[1]);
  double uR = -2;
  
  cout<<"TL = "<<TL<<endl;
  cout<<"TR = "<<TR<<endl;
  cout<<"uL = "<<uL<<endl;
  cout<<"uR = "<<uR<<endl;

  MatrixXcd Gr(L,L);
  MatrixXcd Ga(L,L);
  MatrixXcd Sigma_L(L,L);
  MatrixXcd Sigma_R(L,L);
  MatrixXcd Gamma_L(L,L);
  MatrixXcd Gamma_R(L,L);
  MatrixXcd Tcc(L,L);
  MatrixXcd W(L,L);
  MatrixXcd tmp(L,L);

  Tcc.setZero ();
  for (int i = 0; i < L-1; ++i)
    Tcc(i,i+1) = -h;
  for (int i = 1; i < L; ++i)
    Tcc(i,i-1) = -h;
  
  Sigma_L.setZero ();
  Sigma_R.setZero ();
  Gamma_L.setZero ();
  Gamma_R.setZero ();

  double I = 0;
  double Ie = 0;
  for (int i = 0; i < N; ++i)
    {
      double w = (i-N/2)*dw;
      W.setIdentity ();
      W *= w;
      Sigma_L(0,0) = sfg_solve (w,1);
      Sigma_R(L-1,L-1) = sfg_solve (w,1); 
      Gamma_L(0,0) = -2*imag (Sigma_L(0,0));
      Gamma_R(L-1,L-1) = -2*imag (Sigma_R(L-1,L-1));
      Ga = W-Tcc-(Sigma_L+Sigma_R);
      Gr = Ga.inverse ();
      Ga = Gr.conjugate ();
      tmp = (Gr*Gamma_L*Ga*Gamma_R);
      T[i] = real(tmp.trace());
      assert (T[i] >= 0);
      I += T[i]*(fermi(w,uL,TL) - fermi(w,uR,TR));
      Ie += w*T[i]*(fermi(w,uL,TL) - fermi(w,uR,TR));
    }
  cout<<"number current = "<<I*dw/(2*PI)<<endl;
  cout<<"energy current = "<<Ie*dw/(2*PI)<<endl;
  cout<<"heat current = "<<(Ie-uL*I)*dw/(2*PI)<<endl;
}
