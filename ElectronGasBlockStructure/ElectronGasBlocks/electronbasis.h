#ifndef ELECTRONBASIS_H
#define ELECTRONBASIS_H

#include <armadillo>
#include <iostream>

#define pi 3.1415

using namespace std;
using namespace arma;

class ElectronBasis
{
public:
    mat States; // State n has quantum numbers (e, nx, ny, nz, ms) at row n.

    int Nstates;
    int Nholes;
    int Nshells;
    int NshellsFilled;
    int Nparticles;

    double rs;
    double kstep;
    double L3, L2, L1;

    vec EpsilonMatrix;

    ElectronBasis();
    ElectronBasis(int Nshells_input, int NshellsFilled_input, double rs_input);

    // Declaring functions in class
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s);
    double ei(int q);
    void SetUpEpsilon();
    double epsilon(int i, int j, int a, int b);
    double ReferenceEnergy();

    double Absolute_Difference2(int a, int b);

    // Kroenecker deltas
    int KDelta_integer(int a, int b);
    int KDelta_array(rowvec a, rowvec b);
    int KDelta_k(int a, int b);
    int KDelta_sum(int a, int b, int c, int d);
    int KDelta_spin(int a, int b);

};

#endif // ELECTRONBASIS_H
