#ifndef CCDBLOCKS2_H
#define CCDBLOCKS2_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>
#include <math.h>

#define pi 3.1415

using namespace std;
using namespace arma;

class CCDBlocks2
{
public:
    CCDBlocks2();
    CCDBlocks2(basis_set BASIS);

    basis_set basis;

    int Nmax; int m; int M;

    int Nparticles;
    int Nholes;
    int Nstates;
    int Nchannels;

    int Index(int P, int Sz);

    //
    mat t; mat T; mat X;

    double Abs(double a, double b);
    double CCD();
    double CorrolationEnergy();
    void update_t();

};

#endif // CCDBLOCKS2_H
