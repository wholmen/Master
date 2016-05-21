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
    X1 = zeros<mat>(0,2);
    X2 = zeros<mat>(0,2);

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

void Block::AddCrossStates(rowvec x1, rowvec x2){

    // Function made to add all cross channels

    if ( ! ( any( X1.col(0) == x1(0)) && any( X1.col(1) == x1(1)) )){
        // Adding state
        int n = X1.n_rows;
        X1.insert_rows(n,1);

        X1(n,0) = x1(0); X1(n,1) = x1(1);
    }
    if ( ! ( any( X2.col(0) == x2(0)) && any( X2.col(1) == x2(1)) )){
        // Adding state
        int n = X2.n_rows;
        X2.insert_rows(n,1);

        X2(n,0) = x2(0); X2(n,1) = x2(1);
    }
}

void Block::FinishBlock(){
    // Number of two-state configurations
    Nh = Holes.n_rows;
    Np = Particles.n_rows;
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

    int Nx = X1.n_rows;


    T = zeros<mat>(Nx,Nx);
    V = zeros<mat>(Nx,Nx);

    for (int x1=0; x1<Nx; x1++){
        for (int x2=0; x2<Nx; x2++){

            int a = X1(x1,1); int b = X1(x2,1); int aa = a-Nholes; int bb = b-Nholes;
            int i = X1(x1,0); int j = X1(x2,0);

            T(x1,x2) = t0( Index(i,aa,j,bb,Nholes,Nparticles,Nholes) );
            V(x2,x1) = basis.TwoBodyOperator(i,a,j,b);

        }
    }
}

void Block::SetUpMatrices_Qa(mat &t0){

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

void Block::SetUpMatrices_Qb(mat &t0){

    int Nx = X1.n_rows;

    T = zeros<mat>(Nx,Nx);
    V = zeros<mat>(Nx,Nx);
    T2 = zeros<mat>(Nx,Nx);

    for (int x1=0; x1<Nx; x1++){
        for (int x2=0; x2<Nx; x2++){

            int a = X1(x1,0); int b = X1(x2,0); int aa = a-Nholes; int bb = b-Nholes;
            int i = X1(x1,1); int j = X1(x2,1);

            T(x1,x2) = t0( Index(aa,i,bb,j,Nparticles,Nholes,Nparticles) );
            V(x2,x1) = basis.TwoBodyOperator(a,i,b,j);
            T2(x1,x2) = t0( Index(aa,i,bb,j,Nparticles,Nholes,Nparticles));
        }
    }

    cout << endl << "T: " << endl;
    T.print();

    //cout << endl << "V: " << endl;
    //V.print();


    /*

    mat X1 = zeros<mat>(Np*Nh,2);
    mat X2 = zeros<mat>(Np*Nh,2);

    mat X = zeros<mat>(0,2); mat Xc = zeros<mat>(0,2);
    int n=0;
    for (int ii=0; ii<Nh; ii++){
        for (int aa=0; aa<Np; aa++){

            int i = Holes(ii,0); int a = Particles(aa,0);

            X.insert_rows(n,1); Xc.insert_rows(n,1);
            X(n,0) = i; X(n,1) = a;
            Xc(n,0) = i; Xc(n,1) = i;
            n++;
        }
    }
    // Cross channels set up
    int Nx = X.n_rows;
    cout << Nx << endl;
    T = zeros<mat>(Nx,Nx);
    V = zeros<mat>(Nx,Nx);
    T2 = zeros<mat>(Nx,Nx);

    for (int x1=0; x1<Nx; x1++){
        for (int x2=0; x2<Nx; x2++){

            int i=X(x1,0); int a=X(x1,1); int aa = a-Nholes;
            int j=X(x2,0); int b=X(x2,1); int bb = b-Nholes;

            T(x1,x2) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
        }
    }
    //T.print();

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int n = A + I*Np;

            X1(n,0) = Particles(A,0);
            X1(n,1) = Holes(I,0);
            X2(n,0) = Particles(A,1);
            X2(n,1) = Holes(I,1);
        }
    }

    T = zeros<mat>(Np,Nh);
    V = zeros<mat>(Nh,Np);
    T2 = zeros<mat>(Np,Nh);

    for (int n=0; n<Np*Nh; n++){

        int I = floor(n/Np);
        int A = n % Np;

        int a = X1(n,0); int i = X1(n,1); int b = X2(n,0); int j = X2(n,1);
        int aa = a-Nholes; int bb = b-Nholes;

        T(A,I) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
        V(I,A) = basis.TwoBodyOperator(i,j,a,b);
        T2(A,I) = t0( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
    }

    //T2.print();



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
    //T2.print(); */
}

void Block::SetUpMatrices_Qc(mat &t0){

    mat K = zeros<mat>(Nh*Np,1);
    mat Ktb1 = zeros<mat>(Nh*Np,3);

    //int n=0;
    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int n = A + I*Np;
            //cout << "n: " << n << " K: " << Holes(I,0) << " Ktb1: " << Particles(A,0) << Particles(A,1) << Holes(I,1) << endl;

            K(n,0) = Holes(I,0);
            Ktb1(n,0) = Particles(A,0);
            Ktb1(n,1) = Particles(A,1);
            Ktb1(n,2) = Holes(I,1);
        }
    }


    T = zeros<mat>(Np,Nh);
    V = zeros<mat>(Nh,Np);
    T2 = zeros<mat>(Np,Nh);

    for (int I=0; I<Nh; I++){
        for (int A=0; A<Np; A++){

            int n = A + I*Np;

            T(A,I) = t0( Index( Ktb1(n,0)-Nholes, Ktb1(n,1)-Nholes, Ktb1(n,2), K(n,0), Nparticles,Nparticles,Nholes ));
            V(I,A) = basis.TwoBodyOperator( K(n,0), Ktb1(n,2), Ktb1(n,0), Ktb1(n,1) );
            T2(A,I) = t0( Index( Ktb1(n,0)-Nholes, Ktb1(n,1)-Nholes, Ktb1(n,2), K(n,0), Nparticles,Nparticles,Nholes ));
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
















