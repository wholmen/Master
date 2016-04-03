#ifndef NAIVESOLVERS_H
#define NAIVESOLVERS_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>

#define pi 3.1415

using namespace std;
using namespace arma;

class NaiveSolvers
{
private:
    basis_set basis;
    int Nparticles;
    int Nholes;
    int Nstates;

public:
    NaiveSolvers();
    NaiveSolvers(basis_set BASIS);

    // ------------------
    // Coupled Cluster
    // ------------------

    // Variables
    mat t0, t;

    // Functions
    double CCD();
    double dE_CCD();
    void update_t0();
    void update_t(int degree);


    // Many-Body Perturbation Theory
    double MBPT2_MHJ();

    double MBPT2();
    double MBPT3();
    double MBPT4();
};

#endif // NAIVESOLVERS_H
