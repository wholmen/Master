#include "block.h"

Block::Block(){}

Block::Block(ElectronBasis BASIS, int NHOLES, int NPARTICLES){
    // Block has been initialized
    basis = BASIS;
    Nholes = NHOLES;
    Nparticles = NPARTICLES;

    Holes = zeros<mat>(0,3); // Rows will be added when needed.
    Particles = zeros<mat>(0,3);
    Xph = zeros<mat>(0,3);
    Xhp = zeros<mat>(0,3);
    Kh = zeros<mat>(0,2);
    Khpp = zeros<mat>(0,4);
}

void Block::AddStates(rowvec I, rowvec A){
    // The following if-test is not a beauty.
    // Explained in short terms: It checks wether I have added the samme combination of i,j before. i.e. is the two-hole state already added?
    // I will only need one of each two-hole state
    if ( ! (any( Holes.col(0) == I(0) ) && any( Holes.col(1) == I(1)))  ) { // We have not already stored this exact two-hole state

        int n = Holes.n_rows;
        Holes.insert_rows(n,1);

        Holes(n,0) = I(0); Holes(n,1) = I(1); Holes(n,2) = I(2);
    }
    if ( ! (any( Particles.col(0) == A(0) ) && any( Particles.col(1) == A(1)))  ) { // We have not already stored this exact two-hole state

        int n = Particles.n_rows;
        Particles.insert_rows(n,1);

        Particles(n,0) = A(0); Particles(n,1) = A(1); Particles(n,2) = A(2);
    }
}

void Block::AddPHHPStates(rowvec x1, rowvec x2){
    if ( ! (any( Xph.col(0) == x1(0)) && any( Xph.col(1) == x1(1)) )){

        int n = Xph.n_rows;
        Xph.insert_rows(n,1);

        Xph(n,0) = x1(0); Xph(n,1) = x1(1); Xph(n,2) = x1(2);
    }
    if ( ! (any( Xhp.col(0) == x2(0)) && any( Xhp.col(1) == x2(1)) )){

        int n = Xhp.n_rows;
        Xhp.insert_rows(n,1);

        Xhp(n,0) = x2(0); Xhp(n,1) = x2(1); Xhp(n,2) = x2(2);
    }
}

void Block::AddTripleStates(rowvec k1, rowvec k2){
    //k2.print();
    if ( ! (any( Kh.col(0) == k1(0)))){

        int n = Kh.n_rows;
        Kh.insert_rows(n,1);
        Kh(0) = k1(0); Kh(1) = k1(1);
    }
    if ( ! (any( Khpp.col(0) == k2(0)) && any( Khpp.col(1) == k2(1)) && any( Khpp.col(2) == k2(2)) )){

        int n = Khpp.n_rows;
        Khpp.insert_rows(n,1);
        Khpp(n,0) = k2(0); Khpp(n,1) = k2(1); Khpp(n,2) = k2(2); Khpp(n,3) = k2(3);
    }
}

void Block::FinishBlock(){
    // Number of two-state configurations
    Nh = Holes.n_rows;
    Np = Particles.n_rows;
    Nph = Xph.n_rows;
    Nkh = Kh.n_rows;
    Nkhpp = Khpp.n_rows;
}


void Block::SetUpMatrices_Energy(mat &t){
    // It is important that I do not copy the vector t0 into this function, as that might be disastrous for the memory.
    // I will only send in the memory adress.

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
    V = zeros(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            V(A,I) = basis.TwoBodyOperator( Particles(A,0), Particles(A,1), Holes(I,0), Holes(I,1) );
        }
    }
}

void Block::SetUpMatrices_La(mat &t0){

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

void Block::SetUpMatrices_Lc(mat &t0){

    V = zeros<mat>(Nph,Nph);
    T = zeros<mat>(Nph,Nph);

    for (int x1=0; x1<Nph; x1++){
        for (int x2=0; x2<Nph; x2++){

            int a=Xph(x1,0); int i=Xph(x1,1);
            int b=Xph(x2,0); int j=Xph(x2,1);

            V(x2,x1) = basis.TwoBodyOperator(j,a,b,i);
        }
        for (int x2=0; x2<Nph; x2++){

            int i=Xph(x1,1); int j=Xhp(x2,0);
            int a=Xph(x1,0); int b=Xhp(x2,1);
            int aa=a-Nholes; int bb=b-Nholes;

            T(x1,x2) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
        }
    }
}

void Block::SetUpMatrices_Qa(mat &t0){

    T = zeros(Np,Nh);
    V = zeros(Nh,Np);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int i = Holes(I,0); int j = Holes(I,1);
            int a = Particles(A,0); int b = Particles(A,1);
            int aa = a - Nholes; int bb = b - Nholes;

            T(A,I) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) );
            V(I,A) = basis.TwoBodyOperator(i,j,a,b);
        }
    }
}

void Block::SetUpMatrices_Qb(mat &t0){

    T = zeros<mat>(Nph,Nph);
    V = zeros<mat>(Nph,Nph);

    for (int x1=0; x1<Nph; x1++){
        for (int x2=0; x2<Nph; x2++){

            int i=Xph(x1,1); int j=Xhp(x2,0);
            int a=Xph(x1,0); int b=Xhp(x2,1);
            int aa=a-Nholes; int bb=b-Nholes;

            T(x1,x2) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
            V(x2,x1) = basis.TwoBodyOperator(j,i,b,a);
        }
    }
}

void Block::SetUpMatrices_Qc(mat &t0){

    T = zeros<mat>(Nkhpp, Nkh);
    V = zeros<mat>(Nkh, Nkhpp);

    for (int k1=0; k1<Nkh; k1++){
        for (int k2=0; k2<Nkhpp; k2++){

            int i=Kh(k1,0);   int j=Khpp(k2,0);
            int a=Khpp(k2,1); int b=Khpp(k2,2);
            int aa=a-Nholes;  int bb=b-Nholes;

            T(k2,k1) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
            V(k1,k2) = basis.TwoBodyOperator(i,j,a,b);
        }
    }
}

void Block::SetUpMatrices_Qd(mat &t0){

    mat K = zeros<mat>(Nh*Np,1);
    mat Ktb1 = zeros<mat>(Nh*Np,3);

    //int n=0;
    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int n = A + I*Np;
            //cout << "n: " << n << " K: " << Holes(I,0) << " Ktb1: " << Particles(A,0) << Particles(A,1) << Holes(I,1) << endl;

            K(n,0) = Particles(A,1);
            Ktb1(n,0) = Holes(I,0);
            Ktb1(n,1) = Holes(I,1);
            Ktb1(n,2) = Particles(A,0);
        }
    }


    T = zeros<mat>(Np,Nh);
    V = zeros<mat>(Nh,Np);
    T2 = zeros<mat>(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int n = A + I*Np;

            T(A,I) = t0( Index( Ktb1(n,0), Ktb1(n,1), Ktb1(n,2)-Nholes, K(n,0)-Nholes, Nholes,Nholes,Nparticles ));
            V(I,A) = basis.TwoBodyOperator( K(n,0), Ktb1(n,2), Ktb1(n,0), Ktb1(n,1) );
            T2(A,I) = t0( Index( Ktb1(n,0), Ktb1(n,1), Ktb1(n,2)-Nholes, K(n,0)-Nholes, Nholes,Nholes,Nparticles ));
        }
    }
}

int Block::Index(int p, int q, int r, int s, int Np, int Nq, int Nr){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}
















