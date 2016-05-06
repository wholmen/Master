#include "ccdblocks2.h"

CCDBlocks2::CCDBlocks2()
{

}

CCDBlocks2::CCDBlocks2(basis_set BASIS){

    basis = BASIS;
    Nholes = basis.Nparticles;
    Nholes2 = Nholes*Nholes;
    Nholes3 = Nholes2*Nholes;

    Nstates = basis.nstates;

    Nparticles = Nstates - Nholes;
    Nparticles2 = Nparticles*Nparticles;
    Nparticles3 = Nparticles2*Nparticles;

    //t0 = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);
    //t = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);
    t0 = zeros<vec>(Nparticles2 * Nholes2);
    t = zeros<vec>(Nparticles2 * Nholes2);

}

void CCDBlocks2::Update_Channels(){
    T = zeros<mat>(0,3); int n=0;

    for (int p=0; p<Nstates; p++){
        for (int q=0; q<Nstates; q++){

            if (p != q) { // Pauli exclusion satisfied

                int Nx = basis.states(p,1) + basis.states(q,1);
                int Ny = basis.states(p,2) + basis.states(q,2);
                int Nz = basis.states(p,3) + basis.states(q,3);

                int Sz = basis.states(p,4) + basis.states(p,4);

                int Nmax = basis.Nshells;
                m = 2 * floor(sqrt(Nmax));
                M = 2*m + 1;

                T.insert_rows(n,1);

                T(n,0) = p; T(n,1) = q; T(n,2) = Index(Nx,Ny,Nz,Sz);

                n++;
            }
        }
    }
    Nchannels = T.n_rows;
}

int CCDBlocks2::Index(int Nx, int Ny, int Nz, int Sz){
    return 2*(Nx + m)*M*M*M + 2*(Ny + m)*M*M + 2*(Nz + m)*M + 2*(Sz + 1);
}


double CCDBlocks2::Abs(double a, double b){
    return sqrt( pow(a-b, 2) );
}

double CCDBlocks2::CCD(){
    int n = 0;
    double error = 0.00001;

    Update_Channels();
    update_t();

    double Eold = 0.0;
    double Enew = CorrolationEnergy();

    while ( n < 10 && Abs(Enew,Eold) > error){

        update_t();
        Eold = Enew;
        Enew = CorrolationEnergy();

        n ++;
    }
    cout << "Coupled Cluster method used " << n << " iterations." << endl;
    return Enew;
}





double CCDBlocks2::CorrolationEnergy_NAIVE(){
    // FIRST ATTEMPT AT USING BLOCK METHOD FOR COMPUTING CORROLATION ENERGY

    // NOTE: This function only works for Ns = 2, Nh = 2;
    // This is because no general method is used. Sizes of matrices are hardcoded for
    // self educating purposes.

    // The method is working, meaning that it produces the same results as CorrolationEnergy()
    // This will be the naive approach giving inspiration to a general blocking corrolation energy
    // working for all variations on particles and holes.

    double E = 0.0; // Corrolation energy

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Set up two-particle configuration matrix
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mat Holes = zeros<mat>(0,3);  // A matrix for all the two-particle configurations of holes
    mat Particles = zeros<mat>(0,3); // A matrix for all the two-particle configurations of particles

    // Setting up variables to calculate unique identifier for block computing
    int Nmax = basis.Nshells; m = 2 * floor(sqrt(Nmax)); M = 2*m + 1;


    int n=0; // n will count how many two-state combinations we find. Used as indice in the matrix
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            if (i != j){ // Pauli principle demands that the particles must be unequal

                int Nx = basis.states(i,1) + basis.states(j,1); // Combining x-momentum
                int Ny = basis.states(i,2) + basis.states(j,2); // Combining y-momentum
                int Nz = basis.states(i,3) + basis.states(j,3); // Combining z-momentum
                int Sz = basis.states(i,4) + basis.states(j,4); // Combining spin

                // Adding a new two-hole-state configuration to matrix. (i, j, Identifier)
                Holes.insert_rows(n,1);
                Holes(n,0) = i; Holes(n,1) = j; Holes(n,2) = Index(Nx,Ny,Nz,Sz);
                n++;
            }
        }
    }

    // Setting up particle matrix by the same principles as hole matrix.
    n=0;
    for (int aa=0; aa<Nparticles; aa++){
        for (int bb=0; bb<Nparticles; bb++){
            int a=aa+Nholes; int b=bb+Nholes;

            if (aa != bb){

                int Nx = basis.states(a,1) + basis.states(b,1);
                int Ny = basis.states(a,2) + basis.states(b,2);
                int Nz = basis.states(a,3) + basis.states(b,3);
                int Sz = basis.states(a,4) + basis.states(b,4);

                Particles.insert_rows(n,1);

                Particles(n,0) = a; Particles(n,1) = b; Particles(n,2) = Index(Nx,Ny,Nz,Sz);
                n++;
            }
        }
    }
    int NPARTICLES = Particles.n_rows; int NHOLES = Holes.n_rows; // Number of two-particle configurations


    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Calculate energy
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // For Nh=2; Ns = 2; we will only have 1 block. i+j, a+b will only have one common identifier.
    // As a result, a naive approach can be used by setting up a single matrix for both the interaction and amplitudes
    // I have counted the number of two-particle configurations with the correct identifier to be 12

    mat v = zeros<mat>(2,12); int nv0 = 0; int nv1 = 0;
    mat T = zeros<mat>(12,2); // Matrices for interaction and amplitudes

    for (int i=0; i<NHOLES; i++){
        for (int a=0; a<NPARTICLES; a++){

            if (Holes(i,2) == Particles(a,2)){ // Check if they belong to the same block. i.e. conservation requirements are fulfilled.

                // Filling up the matrices with the correct values. Note again, this is only valid for Ns=2; Nh=2
                if (i == 0){
                    v(0,nv0) = basis.TwoBodyOperator(Particles(a,0),Particles(a,1), Holes(i,0),Holes(i,1));

                    T(nv0,0) = t0( t_index( Particles(a,0)-Nholes,Particles(a,1)-Nholes,Holes(i,0),Holes(i,1) ) );

                    nv0 ++;
                }
                if (i == 1){
                    v(1,nv1) = basis.TwoBodyOperator(Holes(i,0),Holes(i,1),Particles(a,0),Particles(a,1));

                    T(nv1,1) = t0( t_index( Particles(a,0)-Nholes,Particles(a,1)-Nholes,Holes(i,0),Holes(i,1) ) );

                    nv1 ++;
                }
            }
        }
    }
    mat Energy = v*T; // Calculating the energy as a matrix-matrix multiplication

    vec eigval = eig_sym(Energy);

    E = eigval.min();

    return E / 4.0;
}


double CCDBlocks2::CorrolationEnergy2(){
    double E = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += basis.TwoBodyOperator(a,b,i,j) * t0( t_index( aa,bb,i,j) );
                }
            }
        }
    }
    return E / 4.0;
}



int CCDBlocks2::t_index( int a, int b, int i, int j){
    return a + b*Nparticles + i*Nparticles2 + j*Nparticles2*Nholes;
}


void CCDBlocks2::update_t(){
    t0 = t;

}
