#ifndef SOLVER_HPP
#define SOLVER_HPP
#define _USE_MATH_DEFINES
#include <cmath>
#include "armadillo"
#include <string>

using namespace arma;
using namespace std;

class Finance
{
    private:
    // Attributes
    vec m_avec;
    int m_agents;
    double m_m_0;
    int m_mcs;
    double m_beta;
    double m_tax;
    double m_tax_or_no;
    double m_variance;
    double m_savings;
    double m_norm;
    double m_alpha;
    int m_counter;
    double m_gamma;
    string m_filename;
    mt19937_64 generator;
    uniform_real_distribution<double> dis;
    mat m_cij;
    02fb96a34e69af8d8c19226fb51c0e2a557149a8

    public:
    //Values under are chosen to be public for testing purposes.

    void Initialize(int mcs, int agents, double m_0, string filename,
                    double tax_or_no, double tax, double savings,
                    double alpha, double gamma);
    //Functions
    void print_omega(string filename);
    void Metropolis();
    double ss(int i, double delta_m);
    void MonteCarlo();
    void print_avg_dist();
    void print_vec(string filename);
    double ran1();
    double p_dist(int i, int j, double m_cij_val);
    void calc_avg_dist();
};

#endif // ISING_MCINT
