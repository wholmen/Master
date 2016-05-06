#include "solver.h"

Solver::Solver(basis_set BASIS)
{
    basis = BASIS;

    // Calculating important variables

    Nholes = basis.Nparticles; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = basis.nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

    int Nmax = basis.Nshells;
    m = 2 * floor(sqrt(Nmax));
    M = 2*m + 1;


    // Setting up two-state configurations

    TwoBodyConfigurations();
    NPARTICLES = Particles.n_rows;
    NHOLES = Holes.n_rows;
    NSTATES = States.n_rows;


    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);
}


double Solver::CCD(){

    // Set up the first calculation for all amplitudes equal to zero
    double E0 = CorrolationEnergy(); // Can be hardcoded to 0 to save computation cost

    // Generate first set of new amplitudes and do second calculation
    UpdateAmplitudes();
    double E1 = CorrolationEnergy();
    cout << "Energy using blocking: " << CorrolationEnergy() << endl;
    cout << "Energy using naive approach: " << CorrolationEnergy2() << endl;

    // Start the iteration process
    NIterations = 0; tolerance = 1e-6;

    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < 10){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrolationEnergy();
        cout << "Energy using blocking: " << CorrolationEnergy() << endl;
        cout << "Energy using naive approach: " << CorrolationEnergy2() << endl;

        NIterations ++;
    }

    return E1;
}


double Solver::CorrolationEnergy(){

    // Counting and storing all unique identifiers that exist for both Holes and Particles
    vec identifiers = zeros<vec>(0); int NId = 0;
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) { // They have matching identifier and belong to the same block.

                bool IDExist = any( identifiers == Holes(I,2) ); // Checking wether the identifier has already been counted

                if ( ! IDExist ){ // We have found a new identifier. Storing it in the list of identifiers

                    identifiers.insert_rows(NId, 1);
                    identifiers(NId) = Holes(I,2);
                    NId ++;
                }
            }
        }
    }

    Block **blocks = new Block*[NId];

    // Initiate all the blocks needed
    for (int i=0; i<NId; i++){
        blocks[i] = new Block(basis, Nholes, Nparticles);
    }


    // Add all two-state configurations to the corresponding blocks
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) {

                uvec indices = find( identifiers == Holes(I,2) );
                int indice = indices(0);

                blocks[indice] ->AddStates( Holes.row(I), Particles.row(A) );
            }
        }
    }

    // Set up matrices for all blocks
    for (int i=0; i<NId; i++){
        blocks[i] ->SetUpMatrices(t0);
    }

    //blocks[0]->V.print();
    //blocks[0]->T.print();

    // Do the matrix-matrix multiplications for all blocks
    double E = 0;
    for (int i=0; i<NId; i++){

        mat Energy = blocks[i]->V * blocks[i]->T;

        vec eigval = eig_sym(Energy);
        E += eigval.min();
    }

    return E/4;
}


double Solver::CorrolationEnergy2(){
    double E = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += basis.TwoBodyOperator(a,b,i,j) * t0( Index( aa,bb,i,j) );
                }
            }
        }
    }
    return E / 4.0;
}



void Solver::UpdateAmplitudes(){

    t0 = t; // We need to save the previous amplitudes in t0. They will be used to compute the new amplitudes stored in t

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    t0( Index(aa,bb,i,j) ) = basis.TwoBodyOperator(a,b,i,j) / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
    t = t0;

}



double Solver::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int Solver::Index(int p, int q, int r, int s){
    // Tested: Working
    // For consistency: (a,b,i,j) provide different result from (i,j,a,b).

    int Np, Nq, Nr;

    Np = (p < Nholes) ? Nholes : Nparticles;
    Nq = (q < Nholes) ? Nholes : Nparticles;
    Nr = (r < Nholes) ? Nholes : Nparticles;

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}

void Solver::TwoBodyConfigurations(){

    // Function tested. Function working.

    Holes = zeros<mat>(0,3);
    Particles = zeros<mat>(0,3);
    States = zeros<mat>(0,3);

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
                Holes(n,0) = i; Holes(n,1) = j; Holes(n,2) = Identifier(Nx,Ny,Nz,Sz);
                n++;
            }
        }
    }
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

                Particles(n,0) = a; Particles(n,1) = b; Particles(n,2) = Identifier(Nx,Ny,Nz,Sz);
                n++;
            }
        }
    }
    n=0;
    for (int p=0; p<Nstates; p++){
        for (int q=0; q<Nstates; q++){

            if (p != q) { // Pauli exclusion satisfied

                int Nx = basis.states(p,1) + basis.states(q,1);
                int Ny = basis.states(p,2) + basis.states(q,2);
                int Nz = basis.states(p,3) + basis.states(q,3);
                int Sz = basis.states(p,4) + basis.states(p,4);

                States.insert_rows(n,1);

                States(n,0) = p; States(n,1) = q; States(n,2) = Identifier(Nx,Ny,Nz,Sz);
                n++;
            }
        }
    }
}

double Solver::Identifier(int Nx, int Ny, int Nz, int Sz){
    return 2*(Nx + m)*M*M*M + 2*(Ny + m)*M*M + 2*(Nz + m)*M + 2*(Sz + 1);
}




























