#ifndef SOLVER_H
#define SOLVER_H

#include <armadillo>
#include <iostream>
#include <basis_set.h>
#include <math.h>
#include <block.h>

#define pi 3.1415

using namespace std;
using namespace arma;


class Solver
{
public:
    Solver(basis_set BASIS); // Function to initiate class. Sets up two-body configurations and important variables

    basis_set basis; // Variable to hold information on the basis set


    // Declaring important variables

    int Nparticles, Nparticles2, Nparticles3; // Number of particle states to first, second and third order
    int Nholes, Nholes2, Nholes3; // Number of hole states to first, second and third order
    int Nstates; // Total number of states, Nsates = Nholes + Nparticles;

    int NPARTICLES; // Number of two-particle cofigurations
    int NHOLES; // Number of two-hole configurations
    int NSTATES; // Number of two-state configurations
    int NX; // Number of cross states

    int Nmax, m, M; // Variables to create unique identifier for two-body configuration

    double weight; // Weighting the newly calculated amplitudes
    double tolerance; // Convergence criteria for iteration process for Amplitudes
    int NIterations;  // Number of iterations used for calculating amplitudes

    // Declaring matrices
    vec t0; // Amplitudes for iteration n
    vec t ; // Amplitudes for iteration n+1
    mat Holes; // Matrix for all two-hole configurations. Each row contains (state1, state2, identifier)
    mat Particles; // Matrix for all two-particle configurations. Each row contains (state1, state2, identifier)
    mat States; // Matrix for all two-state configurations. Each row contains (state1, state2, identifier)
    mat Xph, Xhp;

    // Diagrams
    void DiagramL0();

    void DiagramLa();
    void DiagramLb();
    void DiagramLc();

    void DiagramQa();
    void DiagramQb();
    void DiagramQc();
    void DiagramQd();

    // Block list
    void CreateBlocks();
    Block **blockspphh; void CreateBlocksPPHH(); int Npphh;
    Block **blocksphph; void CreateBlocksPHPH(); int Nphph;
    Block **blockshphp; void CreateBlocksHPHP(); int Nhphp;

    // Declaring functions for program flow

    double CCD(bool naive); // This is the main function in solver. Calling this function will fully compute the corrolation energy using CCD
    double CorrolationEnergy(); // Function used to calculate energy for every iteration in the amplitudes

    void UpdateAmplitudes(); // Function that updates the amplitudes
    void UpdateAmplitudes_Naive();

    void TwoBodyConfigurations(); // A function call that set up the matrices Holes, Particles and States with sizes (NHOLES,3) (NPARTICLES,3) (NSTATES,3). Each row contains (state1, state2, identifier)


    // Various assisting functions

    double AbsoluteDifference(double a, double b); // Returns sqrt( (a-b)^2 )
    double Identifier(int Nx, int Ny, int Nz, int Sz); // Returns the unique identifier for two-state configuration. Returns 2*(Nx + m)*M^3 + 2*(Ny + m)*M^2 + 2*(Nz + m)*M + 2*(Sz + 1);
    int Index(int p, int q, int r, int s, int Np, int Nq, int Nr); // Returns the index p + q*Np + r*Np*Nq + s*Np*Nq*Nr
};

#endif // SOLVER_H
