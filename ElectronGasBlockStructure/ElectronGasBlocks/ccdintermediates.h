#ifndef CCDINTERMEDIATES_H
#define CCDINTERMEDIATES_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>
#include <math.h>

#define pi 3.1415

using namespace std;
using namespace arma;

class CCDIntermediates
{
public:
    CCDIntermediates();
    CCDIntermediates(basis_set BASIS); // Function to initiate class. Sets up two-body configurations and important variables

    basis_set basis; // Variable to hold information on the basis set


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


    // Declaring functions for program flow

    double CCD(int MaxIterations); // This is the main function in solver. Calling this function will fully compute the corrolation energy using CCD
    double CorrelationEnergy(); // Function used to calculate energy for every iteration in the amplitudes
    void UpdateAmplitudes(); // Function that updates the amplitudes
    vec CCD_ReturnAllIterations();

    // Intermediates
    mat I1, I2, I3, I4;
    void UpdateI1();
    void UpdateI2();
    void UpdateI3();
    void UpdateI4();

    // Various assisting functions

    double AbsoluteDifference(double a, double b); // Returns sqrt( (a-b)^2 )
    int Index(int p, int q, int r, int s, int Np, int Nq, int Nr); // Returns the index p + q*Np + r*Np*Nq + s*Np*Nq*Nr
};




#endif // CCDINTERMEDIATES_H
