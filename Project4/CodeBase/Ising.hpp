#ifndef SOLVER_HPP
#define SOLVER_HPP
#define _USE_MATH_DEFINES
#include <cmath>
#include "armadillo"
#include <string>

using namespace arma;
using namespace std;

class Ising
{
    private:
    // Attributes
    double m_M;
    double m_E;
    int m_spins;
    int m_tot_spins;
    int m_cycles;
    int m_mcs;
    vec m_w;
    int m_size;
    vec m_average;
    double m_init_temp;
    double m_init_temp_sq;
    int m_counter;
    string m_filename;
    mt19937_64 generator;
    uniform_real_distribution<double> dis;
    vec m_E_vals;
    double kb = 1.380649e-23;




    public:
    //Values under are chosen to be public for testing purposes.
    mat m_smatrix;
    double m_variance;
    double m_Etotal_average;
    double m_Mabstotal_average;

    void Initialize(int n_spins, int mcs, double init_temp, int param_1);
    //Functions
    void Metropolis();
    void MonteCarloV1();
    void MonteCarloV2(string filename);
    double up_down(double a);
    int periodic(int i, int limit, int add);
    void calc_variance(int stabile_indx);
    double ran1();
    void init_output(string filename);
    void output();
    void print_E_av(int stabile_indx,string filename);
    void tcoutput(string filename);



};

#endif // ISING_MCINT
