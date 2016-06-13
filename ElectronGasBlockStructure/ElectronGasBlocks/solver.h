#ifndef SOLVER_H
#define SOLVER_H

#include <armadillo>
#include <iostream>
#include <electronbasis.h>
#include <math.h>
#include <block.h>

#define pi 3.1415

using namespace std;
using namespace arma;


class Solver
{
public:
    Solver(ElectronBasis BASIS); // Function to initiate class. Sets up two-body configurations and important variables

    ElectronBasis basis; // Variable to hold information on the basis set


    // Declaring important variables

    int Nparticles, Nparticles2, Nparticles3; // Number of particle states to first, second and third order
    int Nholes, Nholes2, Nholes3; // Number of hole states to first, second and third order
    int Nstates; // Total number of states, Nsates = Nholes + Nparticles;

    int NPARTICLES; // Number of two-particle cofigurations
    int NHOLES; // Number of two-hole configurations
    int NSTATES; // Number of two-state configurations
    int NX; // Number of cross states
    int NKh; int NKh3; int NKp; int NKp3; // Number of triple-states. Either single hole- or particle-state, or combination of h-p-p or p-h-h

    int Nmax, m, M; // Variables to create unique identifier for two-body configuration

    double weight; // Weighting the newly calculated amplitudes
    double tolerance; // Convergence criteria for iteration process for Amplitudes
    int NIterations;  // Number of iterations used for calculating amplitudes

    // Declaring matrices
    vec t0; // Amplitudes for iteration n
    vec t ; // Amplitudes for iteration n+1
    mat Holes; // Matrix for all two-hole configurations. Each row contains (state1, state2, identifier)
    mat Particles; // Matrix for all two-particle configurations. Each row contains (state1, state2, identifier)
    mat Xph, Xhp; // Cross states. p-h config or h-p config
    mat Kh, Khpp; // Triple state channels. h and p+p-h config
    mat Kp, Kphh; // Triple state channels. p and h+h-p config

    // Diagrams
    void DiagramL0();
    void DiagramLa();

    // Without intermediates
    void DiagramLb();
    void DiagramLc();

    void DiagramQa();
    void DiagramQb();
    void DiagramQc();
    void DiagramQd();

    // With intermediates
    void DiagramI1();
    void DiagramI2();
    void DiagramI3();
    void DiagramI4();

    // Block list
    void CreateBlocks();
    Block **blockspphh; void CreateBlocksPPHH(); int Npphh;
    Block **blocksphhp; void CreateBlocksPHHP(); int Nphhp;
    Block **blockskhpp; void CreateBlocksKHPP(); int Nkhpp;
    Block **blockskphh; void CreateBlocksKPHH(); int Nkphh;


    // Declaring functions for program flow
    double CCD(int MaxIterations); // This is the main function in solver. Calling this function will fully compute the corrolation energy using CCD
    double CorrolationEnergy(); // Function used to calculate energy for every iteration in the amplitudes

    void UpdateAmplitudes(); // Function that updates the amplitudes

    void TwoBodyConfigurations(); // A function call that set up the matrices Holes, Particles and States with sizes (NHOLES,3) (NPARTICLES,3) (NSTATES,3). Each row contains (state1, state2, identifier)


    // Various assisting functions

    double AbsoluteDifference(double a, double b); // Returns sqrt( (a-b)^2 )
    double Identifier(int Nx, int Ny, int Nz, int Sz); // Returns the unique identifier for two-state configuration. Returns 2*(Nx + m)*M^3 + 2*(Ny + m)*M^2 + 2*(Nz + m)*M + 2*(Sz + 1);
    int Index(int p, int q, int r, int s); // Returns the index p + q*Np + r*Np*Nq + s*Np*Nq*Nr

    // Diagrams as matrix-matrix multiplication. Included for educational purposes.
    void DiagramMatrixMatrix_Lc();
    void DiagramMatrixMatrix_Qb();
    void DiagramMatrixMatrix_Qc();
    void DiagramMatrixMatrix_Qd();
};

#endif // SOLVER_H
