#ifndef BLOCK_H
#define BLOCK_H

#include <armadillo>
#include <basis_set.h>

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
    Block(basis_set BASIS, int NHOLES, int NPARTICLES);

    basis_set basis; // Using the basis set to compute the interaction when setting up V
    mat Holes; // Storing the Hole states
    mat Particles; // Storing the particle states

    int Nholes, Nparticles;

    mat V; // Interaction matrix
    mat T; // Amplitude matrix

    void AddStates(rowvec I, rowvec A);
    void SetUpMatrices(mat &t0);
    int Index(int p, int q, int r, int s);
};

#endif // BLOCK_H
