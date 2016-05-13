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

    // The following if-test is not a beauty.
    // Explained in short terms: It checks wether I have added the samme combination of i,j before. i.e. is the two-hole state already added?
    // I will only need one of each two-hole state
    if ( ! (any( Holes.col(0) == I(0) ) && any( Holes.col(1) == I(1)))  ) { // We have not already stored this exact two-hole state

        int n = Holes.n_rows;
        Holes.insert_rows(n,1);

        Holes(n,0) = I(0); Holes(n,1) = I(1);

    }
    if ( ! (any( Particles.col(0) == A(0) ) && any( Particles.col(1) == A(1)))  ) { // We have not already stored this exact two-hole state

        int n = Particles.n_rows;
        Particles.insert_rows(n,1);

        Particles(n,0) = A(0); Particles(n,1) = A(1);

    }
    // Function has been tested. It is working as intended.
}


void Block::SetUpMatrices_Energy(mat &t){
    // It is important that I do not copy the vector t0 into this function, as that might be disastrous for the memory.
    // I will only send in the memory adress.

    int Nh = Holes.n_rows;
    int Np = Particles.n_rows;

    V = zeros(Nh,Np);
    T = zeros(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int i = Holes(I,0); int j = Holes(I,1);
            int a = Particles(A,0); int b = Particles(A,1);
            int aa = a - Nholes; int bb = b - Nholes;

            V(I,A) = basis.TwoBodyOperator( i,j,a,b );
            T(A,I) = t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) );
        }
    }
}

void Block::SetUpMatrices_L0(){

    int Nh = Holes.n_rows;
    int Np = Particles.n_rows;

    V = zeros(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            V(A,I) = basis.TwoBodyOperator( Particles(A,0), Particles(A,1), Holes(I,0), Holes(I,1) );
        }
    }
}

void Block::SetUpMatrices_La(mat &t0){

    int Nh = Holes.n_rows;
    int Np = Particles.n_rows;

    V = zeros(Np,Np);
    T = zeros(Np,Nh);

    for (int A=0; A<Np; A++){

        for (int C=0; C<Np; C++){
            V(A,C) = basis.TwoBodyOperator( Particles(A,0), Particles(A,1), Particles(C,0), Particles(C,1) );
        }
        for (int I=0; I<Nh; I++){

            int i = Holes(I,0); int j = Holes(I,1);
            int a = Particles(A,0); int b = Particles(A,1);
            int aa = a - Nholes; int bb = b - Nholes;

            T(A,I) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) );
        }
    }
}

void Block::SetUpMatrices_Lb(mat &t0){

    int Nh = Holes.n_rows;
    int Np = Particles.n_rows;

    V = zeros(Nh,Nh);
    T = zeros(Np,Nh);

    for (int K=0; K<Nh; K++){

        for (int I=0; I<Nh; I++){
            V(K,I) = basis.TwoBodyOperator( Holes(K,0), Holes(K,1), Holes(I,0), Holes(I,1) );
        }
        for (int A=0; A<Np; A++){

            int i = Holes(K,0); int j = Holes(K,1);
            int a = Particles(A,0); int b = Particles(A,1);
            int aa = a - Nholes; int bb = b - Nholes;

            T(A,K) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) );
        }
    }
}

void Block::SetUpMatrices_Qa(mat &t0){

    int Nh = Holes.n_rows;
    int Np = Particles.n_rows;

    T = zeros(Np,Nh);
    V = zeros(Nh,Np);
    T2 = zeros(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int i = Holes(I,0); int j = Holes(I,1);
            int a = Particles(A,0); int b = Particles(A,1);
            int aa = a - Nholes; int bb = b - Nholes;

            T(A,I) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) );

            V(I,A) = basis.TwoBodyOperator(i,j,a,b);

            T2(A,I) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) );
        }
    }
}



int Block::Index(int p, int q, int r, int s, int Np, int Nq, int Nr){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}
















