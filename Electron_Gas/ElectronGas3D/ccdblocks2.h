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

    int Nparticles; int Nparticles2; int Nparticles3;
    int Nholes; int Nholes2; int Nholes3;
    int Nstates;
    int Nchannels;

    vec t0; vec t; mat T; mat X; mat La;

    int Index(int Nx, int Ny, int Nz, int Sz); double m; double M;
    void Update_Channels();

    double Abs(double a, double b);
    double CCD();
    double CorrolationEnergy();
    double CorrolationEnergy_NAIVE();
    double CorrolationEnergy2();
    void update_t();
    void update_t0();
    int t_index(int a, int b, int i, int j);

};

#endif // CCDBLOCKS2_H
