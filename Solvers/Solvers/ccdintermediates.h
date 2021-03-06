#ifndef CCDINTERMEDIATES_H
#define CCDINTERMEDIATES_H

#include <armadillo>
#include <iostream>
#include <../../Pairing_Model/PairingModel/pairingbasis.h>
#include <../../ElectronGasBlockStructure/ElectronGasBlocks/electronbasis.h>
#include <../../NuclearMatter/NuclearMatter/nuclearbasis.h>
#include <math.h>

#define pi 3.14159265359

using namespace std;
using namespace arma;

class CCDIntermediates
{
public:
    CCDIntermediates();
    CCDIntermediates(PairingBasis BASIS); // Function to initiate class. Sets up two-body configurations and important variables
    CCDIntermediates(ElectronBasis BASIS);
    CCDIntermediates(NuclearBasis BASIS);

    PairingBasis pabasis;
    ElectronBasis elbasis;
    NuclearBasis ncbasis;
    int BasisNumber; // A number used to determine wether pabasis or elbasis shall be used inside functions

    // Declaring important variables

    int Nparticles, Nparticles2, Nparticles3; // Number of particle states to first, second and third order
    int Nholes, Nholes2, Nholes3; // Number of hole states to first, second and third order
    int Nstates; // Total number of states, Nsates = Nholes + Nparticles;

    double weight; // Weighting the newly calculated amplitudes
    double tolerance; // Convergence criteria for iteration process for Amplitudes
    int NIterations;  // Number of iterations used for calculating amplitudes


    // Declaring matrices
    vec t0; // Amplitudes for iteration n
    vec t ; // Amplitudes for iteration n+1
    vec EpsilonMatrix;

    // Declaring functions for program flow

    double CCD(int MaxIterations); // This is the main function in solver. Calling this function will fully compute the corrolation energy using CCD
    double CorrelationEnergy(); // Function used to calculate energy for every iteration in the amplitudes
    void UpdateAmplitudes(); // Function that updates the amplitudes
    vec CCD_ReturnAllIterations();
    double v(int p, int q, int r, int s); // Interaction
    double epsilon(int i, int j, int a, int b);
    void SetUpEpsilonMatrix();

    // Intermediates
    mat I1, I2, I3, I4;
    void UpdateI1();
    void UpdateI2();
    void UpdateI3();
    void UpdateI4();

    // Various assisting functions

    double AbsoluteDifference(double a, double b); // Returns sqrt( (a-b)^2 )
    int Index(int p, int q, int r, int s); // Returns the index p + q*Np + r*Np*Nq + s*Np*Nq*Nr
};




#endif // CCDINTERMEDIATES_H
