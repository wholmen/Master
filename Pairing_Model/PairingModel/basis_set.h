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

public:
    mat states;

    int Nshells;
    int Nparticles;
    int nstates;

    double delta;
    double g;

    basis_set();
    basis_set(int nparticles, int nshells, double G, double Delta);

    // Functions for various energy calculations
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s);
    double epsilon(int q);
    double epsilonijab(int i, int j, int a, int b);
    double epsilon4(int i, int j, int k, int l, int a, int b, int c, int d);
    double ReferenceEnergy();

    // Kroenecker Delta's
    int KD_integer(int a, int b);
    int KD_state(int a, int b);
    int KD_spin(int a, int b);

};

#endif // BASIS_SET_H
