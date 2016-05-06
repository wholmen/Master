#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <armadillo>
#include <iostream>

#define pi 3.1415

using namespace std;
using namespace arma;

class basis_set
{
private:
    double rs;
    double kstep;
    double L3, L2, L1;

public:
    mat states; // State n has quantum numbers (e, nx, ny, nz, ms) at row n.

    int nstates;
    int Nparticles;
    int Nshells;

    basis_set();
    basis_set(int nparticles, int nshells, double RS);

    // Declaring functions in class
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s);
    double epsilon(int i, int j, int a, int b);
    double epsilon4(int i, int j, int k, int l, int a, int b, int c, int d);
    double ReferenceEnergy();

    double Absolute_Difference2(int a, int b);

    // Kroenecker deltas
    int KDelta_integer(int a, int b);
    int KDelta_array(rowvec a, rowvec b);
    int KDelta_k(int a, int b);
    int KDelta_sum(int a, int b, int c, int d);
    int KDelta_spin(int a, int b);


//signals:

//public slots:
};

#endif // BASIS_SET_H
