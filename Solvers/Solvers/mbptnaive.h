#ifndef MBPTNAIVE_H
#define MBPTNAIVE_H

#include <armadillo>
#include <iostream>
#include <../../Pairing_Model/PairingModel/pairingbasis.h>
#include <math.h>

#define pi 3.1415

using namespace std;
using namespace arma;


class MBPTNaive
{
public:
    MBPTNaive();
    MBPTNaive(PairingBasis BASIS);

    // Declaring important variables
    PairingBasis basis;
    int Nparticles, Nparticles2, Nparticles3; // Number of particle states to first, second and third order
    int Nholes, Nholes2, Nholes3; // Number of hole states to first, second and third order
    int Nstates; // Total number of states, Nsates = Nholes + Nparticles;

    double MBPT2_MHJ();

    double MBPT2();
    double MBPT3();
    double MBPT4();
};

#endif // MBPTNAIVE_H
