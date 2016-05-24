#ifndef FCI_H
#define FCI_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>
#include <math.h>

#define pi 3.1415

using namespace std;
using namespace arma;


class FCI
{
public:
    FCI();
    FCI(basis_set BASIS);

    basis_set basis;

    int Nparticles, Nparticles2, Nparticles3; // Number of particle states to first, second and third order
    int Nholes, Nholes2, Nholes3; // Number of hole states to first, second and third order
    int Nstates; // Total number of states, Nsates = Nholes + Nparticles;
    double g;

    double CalculateFCI();
    double CalculateCI();
};

#endif // FCI_H
