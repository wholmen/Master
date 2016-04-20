#ifndef CCDBLOCKS_H
#define CCDBLOCKS_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>

#define pi 3.1415

using namespace std;
using namespace arma;

// Creating a class for
class CCDBlocks
{
public:
    CCDBlocks();
    CCDBlocks(basis_set BASIS);

    // Defining variables
    basis_set basis;
    int Nparticles;
    int Nholes;
    int Nstates;

    // Setting up matrices for amplitudes
    mat t, t0;

    // Setting up matrices for diagrams
    mat Vabcd; mat Vabij; mat Vklij; mat Vklcd; mat eps;
    void UpdateVabij();
    void UpdateVabcd();
    void UpdateVklij();
    void UpdateVklcd();

    void UpdateEps();

    // Functions to calculate diagrams


    // Setting up intermediates
    mat I1, I2, I3, I4;
    void UpdateI1();
    void UpdateI2();
    void UpdateI3();
    void UpdateI4();

    // Defining Functions for itarating CCD equations
    double CCD();
    double CorrolationEnergy();
    void update_t();

    // Various functions
    double Abs(double a, double b);
};

#endif // CCDBLOCKS_H
