#include "armadillo"
#include "Black_scholes.hpp"
#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include "time.h"
#include <stdio.h>
#include <tuple>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

// polluting the namespaces
using namespace arma;
using namespace std;

//Initialize the boundary, intial conditions and all parameters needed for calc.
void Black_scholes::Initialize(double T,double L, int N,string filename,
                                double r, double D, double sigma, double E){
  m_filename = filename;
  m_S =vec(N);
  m_utilde=vec(N);
  m_uPrev = vec(N);
  m_x = vec(N);
  m_h = 2*L/((double)N+2);
  m_dt = T/((double)N);
  m_N = N;
  m_T =T;
  m_sigma2 =sigma*sigma/((double)2);
  m_alpha = m_dt/((double)m_h*m_h);
  m_a = (r-D)/((double)sigma*sigma)-1/((double)2);
  m_b = r + m_a*(r-D-0.5*sigma*sigma) - m_a*m_a*sigma*sigma;
  m_E = E;
  m_D = D;
  m_r = r;

  m_x(0)=-L;
  m_x(m_N-1)=L;
  m_utilde(0)=m_uPrev(0)=0.;
  m_S(m_N-1)=E*exp(L);
  m_S(0) = E*exp(-L);
  m_utilde(m_N-1)=m_uPrev(m_N-1)=(m_S(m_N-1)-m_E)*exp(m_a*m_x(m_N-1));

  for(int i= 1;i<m_N-1;i++){
    m_x(i) = -L + (double)i*m_h;
    m_S(i) = E*exp(m_x(i));
    m_uPrev(i)= exp(m_x(i)*m_a)*E*std::max(0.,exp(m_x(i))-1.);
  }
}
//Calculates our V_tilde, aka our know vector in the tridiag algorithm.
void Black_scholes::calc_utilde(double t){
  for(int i=1;i<m_N-1; i++){
    m_utilde(i) = m_sigma2*m_alpha*m_uPrev(i-1) + (2.0-m_sigma2*2.*m_alpha)*m_uPrev(i) +m_sigma2*m_alpha*m_uPrev(i+1);
  }
}
//Used as our main algorithm when running all the functions and impliments the Crank_Nic
//scheme.
void Black_scholes::Crank_Nic(int print_per){
  init_print();
  double t = 0;
  for(int y = 1; y < m_N+1; y++){
    t += m_dt;
    m_utilde(m_N-1)=(m_S(m_N-1)*exp(-t*m_D)-m_E*exp(-m_r*t))*exp(m_a*m_x(m_N-1)+m_b*t);
    calc_utilde(t);
    vec u = Tridiag();
    m_uPrev = u;
    if(y%(m_N/print_per)==0){
      vec V = transform_u_V(u,t);
      print_vals(V,t);
    }
    cout << "\r";
    cout << "Calculated:"<<100*t/m_T <<"%"<<flush;
  }
  cout << " "<<endl;
}
//Impliments the tridiagonal matrix algorithm.
vec Black_scholes::Tridiag(){
  vec d(m_N-1); d.fill(2.+2.*m_alpha*m_sigma2);
  vec b(m_N-2); b.fill(-m_alpha*m_sigma2);
  vec u = m_utilde;
  //Forward eliminate:
  for(int i = 1; i < m_N-2; i++){
    //Normalize row
    b(i-1) /= d(i-1);
    u(i) /= d(i-1);
    d(i-1) = 1.;
    //Eliminate
    u(i+1) += u(i)*m_alpha*m_sigma2;
    d(i) += b(i-1)*m_alpha*m_sigma2;
  }
  b(m_N-3) /=d(m_N-3);
  u(m_N-2) /= d(m_N-3);
  d(m_N-3) = 1.;
  //Backward eliminate
  for(int i = m_N-1; i > 1; i--){
    u(i-1) -= u(i)*b(i-2);
  }
  return u;
}

//Transforms our calculated u(x,\tau) to V(S,t)
vec Black_scholes::transform_u_V(vec u,double t){
  vec V = vec(m_N);
  for(int i =0;i<m_N;i++){
    V(i) = u(i)*exp(-(m_a*m_x(i)+m_b*t));
  }
  return V;
}
//Used to run all funtions to spesifically find the greeks. Could have been made
//smarter, but we decided that it was not a good use of time as it is only
//used to vizuallize the greeks once.
void Black_scholes::Greeks(vec sigma,vec r, string rfilename,string sfilename){
  double r_fast = 0.04;
  double simga_fast = 0.4;
  double T = 1;
  //Note: X and N have to be the same as for main_b
  double X=0.75;
  int N=1e3;
  double D=0.12;
  double E=50;
  double t = 0;
  double t2 = 0;
  int N1 = sigma.n_elem;
  int N2 = r.n_elem;
  string filename="NAN.txt";

  ofstream ofile;
  ofile.open(sfilename);


  for (int j = 0; j < N1; j++){
    ofile << setw(20) << setprecision(8) << sigma(j) << " ";
  }
  ofile << setw(20) << setprecision(8) << " " << endl;

  for (int i = 1; i < 6; i++){
    t = (double)i/5. * T;
    ofile << setw(20) << setprecision(8) << t << " ";
    for(int s =0;s<N1;s++){
      Initialize(T,X,N,filename,r_fast,D,sigma(s),E);
      m_utilde(m_N-1)=(m_S(m_N-1)*exp(-t*m_D)-m_E*exp(-m_r*t))*exp(m_a*m_x(m_N-1)+m_b*t);
      calc_utilde(t);
      vec u = Tridiag();
      vec V = transform_u_V(u,t);
      ofile << setw(20) << setprecision(8) << V(m_N-2) << " ";
    }
    ofile << setw(20) << setprecision(8) <<  " " << endl;
  }
  ofile.close();

  ofile.open(rfilename);
  for (int j = 0; j < N2; j++){
    ofile << setw(20) << setprecision(8) << r(j)<< " ";
  }
  ofile << setw(20) << setprecision(8) << " " << endl;

  for (int i = 1; i < 6; i++){
    t2 = (double)i/5. * T;
    ofile << setw(20) << setprecision(8) << t2 << " ";
    for(int k=0;k<N2;k++){
      Initialize(T,X,N,filename,r(k),D,simga_fast,E);
      m_utilde(m_N-1)=(m_S(m_N-1)*exp(-t2*m_D)-m_E*exp(-m_r*t2))*exp(m_a*m_x(m_N-1)+m_b*t2);
      calc_utilde(t2);
      vec u = Tridiag();
      vec V = transform_u_V(u,t2);

      ofile << setw(20) << setprecision(8) << V(m_N-2) << " ";
    }
    ofile << setw(20) << setprecision(8) <<  " " << endl;
  }
  ofile.close();
}
//Empties the txt file, and adds the S-values at the top.
void Black_scholes::init_print(){
  ofstream ofile;
  ofile.open(m_filename);
  for(int i=0;i<m_N-1;i++){
  ofile << setw(20) << setprecision(8) << m_S(i) << " ";
  }
  ofile << setw(20) << setprecision(8) << m_S(m_N-1)<<endl;
  vec V_0 = transform_u_V(m_uPrev,0);
  print_vals(V_0,0);
  ofile.close();
}
//Print all V values.
void Black_scholes::print_vals(vec u,double t){
  ofstream ofile;
  ofile.open(m_filename, fstream::app);
  ofile << setw(20) << setprecision(8) << t << " ";
  for(int i=0;i<m_N-1; i++){
    ofile << setw(20) << setprecision(8) << u(i) << " ";
  }
  ofile << setw(20) << setprecision(8) << u(m_N-1)<<endl;
  ofile.close();
}
