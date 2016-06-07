#ifndef BLOCK_H
#define BLOCK_H

#include <armadillo>
#include <electronbasis.h>

using namespace arma;

class Block
{
public:

    // Block Class will be initialized as a single block corresponding to one identifier.
    // The class will be used for mainly three purposes
    // 1. Store all two-body combinations of hole and particles corresponding to this block
    // 2. Set up an interaction matrix, V, and amplitude matrix T corresponding to this block
    // 3. Returns the matrices to compute energy or new amplitudes

    Block();
    Block(ElectronBasis BASIS, int NHOLES, int NPARTICLES);

    ElectronBasis basis; // Using the basis set to compute the interaction when setting up V
    mat Holes; // Storing the Hole states
    mat Particles; // Storing the particle states
    mat Xph;
    mat Xhp;
    mat Kh, Khpp;

    int Nholes, Nparticles;
    int Nx; // Number of cross states
    int Nph;
    int Nkh; int Nkhpp;

    mat V; // Interaction matrix
    mat T; // Amplitude matrix
    mat T2; // Amplitude matrix nr two for second order diagrams

    void AddStates(rowvec I, rowvec A);
    void AddPHHPStates(rowvec x1, rowvec x2);
    void AddTripleStates(rowvec k1, rowvec k2);

    int Index(int p, int q, int r, int s, int Np, int Nq, int Nr);

    void FinishBlock();
    int Nh; int Np;

    // Set up matrices for diagram calculations
    void SetUpMatrices_Energy(mat &t);
    void SetUpMatrices_L0();
    void SetUpMatrices_La(mat &t0);
    void SetUpMatrices_Lb(mat &t0);
    void SetUpMatrices_Lc(mat &t0);
    void SetUpMatrices_Qa(mat &t0);
    void SetUpMatrices_Qb(mat &t0); mat Qb;
    void SetUpMatrices_Qc(mat &t0);
    void SetUpMatrices_Qd(mat &t0);
};

#endif // BLOCK_H
