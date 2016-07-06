#ifndef NUCLEARBASIS_H
#define NUCLEARBASIS_H

#include <armadillo>
#include <iostream>

#define pi 3.1415

using namespace std;
using namespace arma;

class NuclearBasis
{
public:
    NuclearBasis();
    NuclearBasis(int Nshells_input, int NFilledShells_input, double rs_input, bool OnlyNeutrons);

    mat States;
    int Nshells; int NfilledShells; int Nholes; int Nstates; int Nstates2; int Nstates3; int Nparticles; int Nparticles2; int Nholes2;

    double rs;
    double kstep;
    double L3, L2, L1;

    // Values for computing Minnesota potential
    double v0R, v0T, v0S, kappaR, kappaT, kappaS;

    // Functions for energy calculations
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s); // Minnesota potential as nuclear interaction
    double Minnesota(int p, int q, int r, int s);
    double ei(int q, vec &v);
    double epsilon(int i, int j, int a, int b, vec &v); // Returns e_i + e_j - e_a - e_b
    double ReferenceEnergy();

    int KD_k(int a, int b);
    double Absolute_Difference2(int a, int b);

    // Kroenecker Delta
    int KD_integer(int a, int b);
    int KD_array(rowvec a, rowvec b);
    int KD_sum(int a, int b, int c, int d);
    int KD_spin(int a, int b);
    int KD_isospin(int a, int b);
};




#endif // NUCLEARBASIS_H

