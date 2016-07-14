#ifndef PAIRINGBASIS_H
#define PAIRINGBASIS_H

#include <armadillo>
#include <iostream>

#define pi 3.14159265359

using namespace std;
using namespace arma;

class PairingBasis
{
private:

public:
    mat States;

    int Nshells;
    int NshellsFilled;
    int Nholes;
    int Nstates;
    int Nparticles;

    vec EpsilonMatrix;

    double delta;
    double g;

    PairingBasis();
    PairingBasis(int Nshells_input, int NshellsFilled_input, double g_input, double delta_input);

    // Functions for various energy calculations
    double OneBodyOperator(int p, int q);
    double TwoBodyOperator(int p, int q, int r, int s);
    double ei(int q);
    double epsilon(int i, int j, int a, int b);
    double epsilon4(int i, int j, int k, int l, int a, int b, int c, int d);
    double ReferenceEnergy();
    void SetUpEpsilon();

    // Kroenecker Delta's
    int KD_integer(int a, int b);
    int KD_state(int a, int b);
    int KD_spin(int a, int b);

};

#endif // PAIRINGBASIS_H
