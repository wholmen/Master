#ifndef CCDINTERMEDIATE_H
#define CCDINTERMEDIATE_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>

#define pi 3.1415

using namespace std;
using namespace arma;

class CCDIntermediate
{
private:
    basis_set basis;
    int Nparticles;
    int Nholes;
    int Nstates;

public:
    CCDIntermediate();
    CCDIntermediate(basis_set BASIS);

    // ------------------
    // Coupled Cluster
    // ------------------

    // Variables
    mat t0, t;

    // Functions
    double CCD();
    double dE_CCD();
    void update_t0();
    void update_t();

    // Intermediates
    mat I1, I2, I3, I4;
    void UpdateI1();
    void UpdateI2();
    void UpdateI3();
    void UpdateI4();
};

#endif // CCDINTERMEDIATE_H
