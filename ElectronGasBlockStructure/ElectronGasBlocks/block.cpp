#include "block.h"

Block::Block()
{

}

Block::Block(basis_set BASIS, int NHOLES, int NPARTICLES){
    // Block has been initialized
    basis = BASIS;
    Nholes = NHOLES;
    Nparticles = NPARTICLES;

    Holes = zeros<mat>(0,2); // Rows will be added when needed.
    Particles = zeros<mat>(0,2);

}


void Block::AddStates(rowvec I, rowvec A){
    int n;

    // The following if-test is not a beauty.
    // Explained in short terms: It checks wether I have added the samme combination of i,j before. i.e. is the two-hole state already added?
    // I will only need one of each two-hole state
    if ( ! (any( Holes.col(0) == I(0) ) && any( Holes.col(1) == I(1)))  ) { // We have not already stored this exact two-hole state

        n = Holes.n_rows;
        Holes.insert_rows(n,1);

        Holes(n,0) = I(0); Holes(n,1) = I(1);

    }

    if ( ! (any( Particles.col(0) == A(0) ) && any( Particles.col(1) == A(1)))  ) { // We have not already stored this exact two-hole state

        n = Particles.n_rows;
        Particles.insert_rows(n,1);

        Particles(n,0) = A(0); Particles(n,1) = A(1);

    }
    // Function has been tested. It is working as intended.

}


void Block::SetUpMatrices(mat &t0){
    // It is important that I do not copy the vector t0 into this function, as that might be disastrous for the memory.
    // I will only send in the memory adress.

    int Nh = Holes.n_rows;
    int Np = Particles.n_rows;

    V = zeros(Nh,Np);
    T = zeros(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            V(I,A) = basis.TwoBodyOperator( Holes(I,0), Holes(I,1), Particles(A,0), Particles(A,1) );
            T(A,I) = t0( Index( Particles(A,0)-Nholes, Particles(A,1)-Nholes, Holes(I,0), Holes(I,1) ) );
        }
    }
}


int Block::Index(int p, int q, int r, int s){
    // Tested: Working
    // For consistency: (a,b,i,j) provide different result from (i,j,a,b).

    int Np, Nq, Nr;

    Np = (p < Nholes) ? Nholes : Nparticles;
    Nq = (q < Nholes) ? Nholes : Nparticles;
    Nr = (r < Nholes) ? Nholes : Nparticles;

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}





