#ifndef BASIS_SET2D_H
#define BASIS_SET2D_H

#include <armadillo>
#include <iostream>

#define pi 3.1415

using namespace std;
using namespace arma;

class basis_set2D
{
public:
private:
    int Nparticles;
    int Nshells;
    int MagicNumber;
    int nstates;
    double rs;
    double kstep;
    double L3, L2, L1;

public:
    mat states;

    basis_set2D();
    basis_set2D(int nparticles, int nshells, double RS);

    // Declaring functions in class
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s);
    double Absolute_Difference2(int a, int b);
    double ReferenceEnergy();
    double ReferenceEnergy2();
    int KDelta_integer(int a, int b);
    int KDelta_array(rowvec a, rowvec b);
    int KDelta_k(int a, int b);
    int KDelta_sum(int a, int b, int c, int d);
    int KDelta_spin(int a, int b);

//signals:

//public slots:
};

#endif // BASIS_SET2D_H
