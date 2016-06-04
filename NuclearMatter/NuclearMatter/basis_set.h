#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <armadillo>
#include <iostream>

#define pi 3.1415

using namespace std;
using namespace arma;

class basis_set
{
public:
    basis_set();
    basis_set(int Nshells_input, int NFilledShells_input, double rs_input);

    mat States;
    int Nshells; int NfilledShells; int Nholes; int Nstates;

    double rs;
    double kstep;
    double L3, L2, L1;

    // Values for computing Minnesota potential
    double v0R, v0T, v0S, kappaR, kappaT, kappaS;

    // Functions for energy calculations
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s); // Minnesota potential as nuclear interaction
    double ei(int q);
    double epsilon(int i, int j, int a, int b); // Returns e_i + e_j - e_a - e_b
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




#endif // BASIS_SET_H

