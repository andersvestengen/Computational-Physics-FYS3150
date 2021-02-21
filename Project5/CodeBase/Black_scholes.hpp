#ifndef SOLVER_HPP
#define SOLVER_HPP
#define _USE_MATH_DEFINES
#include <cmath>
#include "armadillo"
#include <string>

using namespace arma;
using namespace std;

class Black_scholes
{
    private:
      double m_h;
      double m_dt;
      int m_N;
      double m_alpha;
      double m_alpha2;
      string m_filename;
      double m_T;
      vec m_S;
      double m_a;
      double m_b;
      vec m_x ;
      double m_sigma2;
      double m_D;
      vec m_utilde;
      vec m_uPrev;
      double m_E;
      double m_r;
    public:
      void Initialize(double T,double L, int N,string filename,
                      double r, double D, double sigma, double E);
      void D1d_explicit();
      void calc_utilde(double t);
      void print_vals(vec u, double t);
      void init_print();
      void Crank_Nic(int print_per);
      void Greeks(vec sigma,vec r, string rfilename,string sfilename);
      vec transform_u_V(vec u,double t);
      vec Tridiag();
};

#endif // ISING_MCINT
