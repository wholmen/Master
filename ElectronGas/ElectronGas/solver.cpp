#include "solver.h"

Solver::Solver(ElectronBasis BASIS)
{
    basis = BASIS;

    // Calculating important variables
    Nholes = basis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = basis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

    int Nmax = basis.Nshells;
    m = 2 * floor(sqrt(Nmax));
    M = 2*m + 1;

    // Weight when adding diagrams to new amplitudes
    weight = 1.0;

    // Setting up two-state configurations
    TwoBodyConfigurations();
    CreateBlocks();

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);
}


double Solver::CCD(int MaxIterations){
    // Set up the first calculation for all amplitudes equal to zero
    double E0 = CorrolationEnergy(); // Can be hardcoded to 0 to save computation cost

    // Generate first set of new amplitudes and do second calculation
    UpdateAmplitudes();
    double E1 = CorrolationEnergy();

    // Start the iteration process
    NIterations = 0; tolerance = 1e-12;
    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < MaxIterations){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrolationEnergy();
        NIterations ++;
    }
    return E1;
}

vec Solver::CCD_ReturnAllIterations(){
    // Variation of CCD() made to store correlation energy for every iteration
    vec energies = zeros<vec>(0);

    double E0 = CorrolationEnergy();
    UpdateAmplitudes();
    double E1 = CorrolationEnergy();

    NIterations = 0; tolerance = 1e-12;
    energies.insert_rows(NIterations,1);
    energies(NIterations) = E1; NIterations++;

    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < 50){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrolationEnergy();
        energies.insert_rows(NIterations,1);
        energies(NIterations) = E1;

        NIterations ++;
    }
    return energies;
}

double Solver::CorrolationEnergy(){

    // Set up matrices for all blocks
    for (int i=0; i<Npphh; i++){
        blockspphh[i] ->SetUpMatrices_Energy(t);
    }

    // Do the matrix-matrix multiplications for all blocks
    double E = 0;
    for (int i=0; i<Npphh; i++){

        mat Energy = blockspphh[i]->V * blockspphh[i]->T;

        vec eigval = eig_sym(Energy);
        E += sum(eigval);
    }
    return E / 4.0;
}

void Solver::UpdateAmplitudes(){
    t0 = t; // We need to save the previous amplitudes in t0. They will be used to compute the new amplitudes stored in t
    t = zeros(Nparticles2*Nholes2);

    // Block approach without intermediates.
    DiagramL0();
    DiagramLa();
    DiagramQc();
    DiagramQd();

    // Block approach with intermediates
    DiagramI1(); // Lb+Qa
    DiagramI2(); // Lc+Qb

    // Adding weight factor
    if (weight != 0) t = weight*t + (1-weight)*t0;
}

void Solver::CreateBlocks(){

    CreateBlocksPPHH();
    CreateBlocksPHHP();
    CreateBlocksKHPP();
    CreateBlocksKPHH();
}

void Solver::CreateBlocksPPHH(){

    vec identifiers = zeros<vec>(0); Npphh = 0;
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if (Holes(I,2) == Particles(A,2) ) { // Block found

                bool IDExist = any( identifiers == Holes(I,2) );

                if ( ! IDExist){
                    identifiers.insert_rows(Npphh,1);
                    identifiers(Npphh) = Holes(I,2);
                    Npphh ++;
                }
            }
        }
    }
    blockspphh = new Block*[Npphh];

    for (int n=0; n<Npphh; n++){
        blockspphh[n] = new Block(basis, Nholes, Nparticles);
    }

    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if (Holes(I,2) == Particles(A,2)){

                uvec indices = find( identifiers == Holes(I,2) );
                int indice = indices(0);

                blockspphh[indice]->AddStates( Holes.row(I), Particles.row(A) );
            }
        }
    }
    for (int n=0; n<Npphh; n++){
        blockspphh[n]->FinishBlock();
    }

    for (int n=0; n<Npphh; n++){
        blockspphh[n]->Epsilonpphh();
    }
}

void Solver::CreateBlocksPHHP(){

    vec identifiers = zeros<vec>(0); Nphhp = 0;
    for (int x1=0; x1<NX; x1++){
        for (int x2=0; x2<NX; x2++){

            if (Xph(x1,2) == Xhp(x2,2)){

                bool IDExist = any( identifiers == Xph(x1,2) );

                if ( ! IDExist){
                    identifiers.insert_rows(Nphhp,1);
                    identifiers(Nphhp) = Xph(x1,2);
                    Nphhp ++;
                }
            }
        }
    }
    blocksphhp = new Block*[Nphhp];

    for (int n=0; n<Nphhp; n++) blocksphhp[n] = new Block(basis,Nholes,Nparticles);

    for (int x1=0; x1<NX; x1++){
        for (int x2=0; x2<NX; x2++){

            if (Xph(x1,2) == Xhp(x2,2)){

                uvec indices = find( identifiers == Xph(x1,2) );
                int indice = indices(0);

                blocksphhp[indice]->AddCrossStates(Xph.row(x1),Xhp.row(x2));
            }
        }
    }
    for (int n=0; n<Nphhp; n++) blocksphhp[n]->FinishBlock();

    for (int n=0; n<Nphhp; n++) blocksphhp[n]->Epsilonphhp();
}

void Solver::CreateBlocksKHPP(){
    vec identifiers = zeros<vec>(0); Nkhpp = 0;
    for (int k1=0; k1<NKh; k1++){
        for (int k2=0; k2<NKh3; k2++){

            if (Kh(k1,1) == Khpp(k2,3) ) { // Block found

                bool IDExist = any( identifiers == Kh(k1,1) );

                if ( ! IDExist){
                    identifiers.insert_rows(Nkhpp,1);
                    identifiers(Nkhpp) = Kh(k1,1);
                    Nkhpp ++;
                }
            }
        }
    }
    blockskhpp = new Block*[Nkhpp];

    for (int n=0; n<Nkhpp; n++) blockskhpp[n] = new Block(basis,Nholes,Nparticles);

    for (int k1=0; k1<NKh; k1++){
        for (int k2=0; k2<NKh3; k2++){

            if (Kh(k1,1) == Khpp(k2,3) ) { // Block found

                uvec indices = find( identifiers == Kh(k1,1));
                int indice = indices(0);

                blockskhpp[indice]->AddTripleStates( Kh.row(k1), Khpp.row(k2) );
            }
        }
    }
    for (int n=0; n<Nkhpp; n++) blockskhpp[n]->FinishBlock();

    for (int n=0; n<Nkhpp; n++) blockskhpp[n]->Epsilonkhpp();
}

void Solver::CreateBlocksKPHH(){
    vec identifiers = zeros<vec>(0); Nkphh = 0;
    for (int k1=0; k1<NKp; k1++){
        for (int k2=0; k2<NKp3; k2++){

            if (Kp(k1,1) == Kphh(k2,3) ) { // Block found

                bool IDExist = any( identifiers == Kp(k1,1) );

                if ( ! IDExist){
                    identifiers.insert_rows(Nkphh,1);
                    identifiers(Nkphh) = Kp(k1,1);
                    Nkphh ++;
                }
            }
        }
    }
    blockskphh = new Block*[Nkphh];

    for (int n=0; n<Nkphh; n++) blockskphh[n] = new Block(basis,Nholes,Nparticles);

    for (int k1=0; k1<NKp; k1++){
        for (int k2=0; k2<NKp3; k2++){

            if (Kp(k1,1) == Kphh(k2,3) ) { // Block found

                uvec indices = find( identifiers == Kp(k1,1));
                int indice = indices(0);

                blockskphh[indice]->AddTripleStates( Kp.row(k1), Kphh.row(k2) );
            }
        }
    }
    for (int n=0; n<Nkphh; n++) blockskphh[n]->FinishBlock();

    for (int n=0; n<Nkphh; n++) blockskphh[n]->Epsilonkphh();
}


void Solver::DiagramL0(){


    // Set up matrices for all blockspphh
    for (int i=0; i<Npphh; i++){
        blockspphh[i] ->SetUpMatrices_L0();
    }
    // For L0 diagram, not computation is needed. Only aligning V elements to t.

    //Align elements to t.
    for (int n=0; n<Npphh; n++){

        mat L0 = blockspphh[n]->V / blockspphh[n]->epsilon;

        for (double I=0; I < blockspphh[n]->Holes.n_rows ; I++){
            for (double A=0; A < blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);

                t( Index(a,b,i,j) ) += L0(A,I);
            }
        }
    }
}

void Solver::DiagramLa(){

    // Tested. Working.

    for (int i=0; i<Npphh; i++){
        blockspphh[i]->SetUpMatrices_La(t0);
    }

    for (int n=0; n<Npphh; n++){
        mat La = 0.5*blockspphh[n]->V * blockspphh[n]->T / blockspphh[n]->epsilon;

        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);

                t( Index(a,b,i,j) ) += La(A,I);
            }
        }
    }
}

void Solver::DiagramQc(){
    // Block Approach
    for (int n=0; n<Nkhpp; n++) blockskhpp[n]->SetUpMatrices_Qc(t0);

    for (int n=0; n<Nkhpp; n++){

        mat Qc = 0.5*blockskhpp[n]->T * blockskhpp[n]->V * blockskhpp[n]->T2 / blockskhpp[n]->epsilon;

        for (int k1=0; k1<blockskhpp[n]->Nk1; k1++){
            for (int k2=0; k2<blockskhpp[n]->Nk3; k2++){

                int j=blockskhpp[n]->K1(k1,0); int i=blockskhpp[n]->K3(k2,0);
                int a=blockskhpp[n]->K3(k2,1); int b=blockskhpp[n]->K3(k2,2);

                t( Index(a,b,i,j)) -= Qc( k2,k1 );
                t( Index(a,b,j,i)) += Qc( k2,k1 );
            }
        }
    }
}

void Solver::DiagramQd(){
    // Block Approach
    for (int n=0; n<Nkphh; n++) blockskphh[n]->SetUpMatrices_Qd(t0);

    for (int n=0; n<Nkphh; n++){

        mat Qd = 0.5*blockskphh[n]->T * blockskphh[n]->V * blockskphh[n]->T2 / blockskphh[n]->epsilon;

        for (int k1=0; k1<blockskphh[n]->Nk1; k1++){
            for (int k2=0; k2<blockskphh[n]->Nk3; k2++){

                int a=blockskphh[n]->K1(k1,0); int b=blockskphh[n]->K3(k2,0);
                int i=blockskphh[n]->K3(k2,1); int j=blockskphh[n]->K3(k2,2);

                t( Index(a,b,i,j)) -= Qd( k2,k1 );
                t( Index(b,a,i,j)) += Qd( k2,k1 );
            }
        }
    }
}

void Solver::DiagramI1(){
    for (int n=0; n<Npphh; n++) blockspphh[n]->SetUpMatrices_I1(t0);

    for (int n=0; n<Npphh; n++){

        mat I1 = 0.5* blockspphh[n]->T * blockspphh[n]->I1 / blockspphh[n]->epsilon;

        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);

                t( Index(a,b,i,j) ) += I1(A,I);
            }
        }
    }
}

void Solver::DiagramI2(){
    for (int n=0; n<Nphhp; n++) blocksphhp[n]->SetUpMatrices_I2(t0);

    for (int n=0; n<Nphhp; n++){

        mat I2 = blocksphhp[n]->T * blocksphhp[n]->I2 / blocksphhp[n]->epsilon;

        for (int x1=0; x1<blocksphhp[n]->Nph; x1++){
            for (int x2=0; x2<blocksphhp[n]->Nph; x2++){

                int i=blocksphhp[n]->Xph(x1,1); int j=blocksphhp[n]->Xhp(x2,0);
                int a=blocksphhp[n]->Xph(x1,0); int b=blocksphhp[n]->Xhp(x2,1);

                t( Index(a,b,i,j) ) += I2(x1,x2);
                t( Index(a,b,j,i) ) -= I2(x1,x2);
                t( Index(b,a,i,j) ) -= I2(x1,x2);
                t( Index(b,a,j,i) ) += I2(x1,x2);
            }
        }
    }
}

// Following functions are various assisting functions for program flow.

double Solver::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int Solver::Index(int a, int b, int i, int j){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state
    return (a-Nholes) + (b-Nholes)*Nparticles + i*Nparticles2 + j*Nparticles2*Nholes;
}

void Solver::TwoBodyConfigurations(){

    // Function tested. Function working.

    Holes = zeros<mat>(0,3);
    Particles = zeros<mat>(0,3);

    int n=0; // n will count how many two-state combinations we find. Used as indice in the matrix
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            if (i != j){ // Pauli principle demands that the particles must be unequal

                // Setting up direct channels for holes

                // Two-hole momentum and spin
                int Nx = basis.States(i,1) + basis.States(j,1); // Combining x-momentum
                int Ny = basis.States(i,2) + basis.States(j,2); // Combining y-momentum
                int Nz = basis.States(i,3) + basis.States(j,3); // Combining z-momentum
                int Sz = basis.States(i,4) + basis.States(j,4); // Combining spin

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

                int Nx = basis.States(a,1) + basis.States(b,1);
                int Ny = basis.States(a,2) + basis.States(b,2);
                int Nz = basis.States(a,3) + basis.States(b,3);
                int Sz = basis.States(a,4) + basis.States(b,4);

                Particles.insert_rows(n,1);
                Particles(n,0) = a; Particles(n,1) = b; Particles(n,2) = Identifier(Nx,Ny,Nz,Sz);
                n++;
            }
        }
    }
    NPARTICLES = Particles.n_rows;
    NHOLES = Holes.n_rows;


    Xhp = zeros<mat>(0,3);
    Xph = zeros<mat>(0,3);

    n=0;
    for (int i=0; i<Nholes; i++){
        for (int aa=0; aa<Nparticles; aa++){
            int a = aa + Nholes;

            // Two-state
            int Nxc = basis.States(i,1) - basis.States(a,1);
            int Nyc = basis.States(i,2) - basis.States(a,2);
            int Nzc = basis.States(i,3) - basis.States(a,3);
            int Szc = basis.States(i,4) - basis.States(a,4);

            Xhp.insert_rows(n,1);
            Xph.insert_rows(n,1);

            Xhp(n,0) = i; Xhp(n,1) = a; Xhp(n,2) = Identifier(Nxc,Nyc,Nzc,Szc);
            Xph(n,0) = a; Xph(n,1) = i; Xph(n,2) = Identifier(-Nxc,-Nyc,-Nzc,-Szc);

            n++;
        }
    }
    NX = Xhp.n_rows;

    Kh = zeros<mat>(0,2); Khpp = zeros<mat>(0,4);

    n=0;
    for (int i=0; i<Nholes; i++){

        int Nx = basis.States(i,1); // Combining x-momentum
        int Ny = basis.States(i,2); // Combining y-momentum
        int Nz = basis.States(i,3); // Combining z-momentum
        int Sz = basis.States(i,4); // Combining spin

        // Adding a new two-hole-state configuration to matrix. (i, j, Identifier)
        Kh.insert_rows(n,1);
        Kh(n,0) = i; Kh(n,1) = Identifier(Nx,Ny,Nz,Sz);
        n++;
    }
    n=0;
    for (int i=0; i<Nholes; i++){
        for (int aa=0; aa<Nparticles; aa++){
            for (int bb=0; bb<Nparticles; bb++){

                int a=aa+Nholes; int b=bb+Nholes;
                if (aa != bb){
                    int Nx = basis.States(a,1) + basis.States(b,1) - basis.States(i,1);
                    int Ny = basis.States(a,2) + basis.States(b,2) - basis.States(i,2);
                    int Nz = basis.States(a,3) + basis.States(b,3) - basis.States(i,3);
                    int Sz = basis.States(a,4) + basis.States(b,4) - basis.States(i,4);

                    Khpp.insert_rows(n,1);
                    Khpp(n,0) = i; Khpp(n,1) = a; Khpp(n,2) = b; Khpp(n,3) = Identifier(Nx,Ny,Nz,Sz);
                }
            }
        }
    }
    NKh3 = Khpp.n_rows; NKh = Kh.n_rows;


    Kp = zeros<mat>(0,2); Kphh = zeros<mat>(0,4);
    n=0;
    for (int aa=0; aa<Nparticles; aa++){
        int a = aa+Nholes;

        int Nx = basis.States(a,1); // Combining x-momentum
        int Ny = basis.States(a,2); // Combining y-momentum
        int Nz = basis.States(a,3); // Combining z-momentum
        int Sz = basis.States(a,4); // Combining spin

        // Adding a new two-hole-state configuration to matrix. (i, j, Identifier)
        Kp.insert_rows(n,1);
        Kp(n,0) = a; Kp(n,1) = Identifier(Nx,Ny,Nz,Sz);
        n++;
    }
    n=0;
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){

                int a=aa+Nholes;
                if (i != j){
                    int Nx = basis.States(i,1) + basis.States(j,1) - basis.States(a,1);
                    int Ny = basis.States(i,2) + basis.States(j,2) - basis.States(a,2);
                    int Nz = basis.States(i,3) + basis.States(j,3) - basis.States(a,3);
                    int Sz = basis.States(i,4) + basis.States(j,4) - basis.States(a,4);

                    Kphh.insert_rows(n,1);
                    Kphh(n,0) = a; Kphh(n,1) = i; Kphh(n,2) = j; Kphh(n,3) = Identifier(Nx,Ny,Nz,Sz);
                    n++;
                }
            }
        }
    }
    NKp3 = Kphh.n_rows; NKp = Kp.n_rows;
}

double Solver::Identifier(int Nx, int Ny, int Nz, int Sz){
    return 2*(Nx + m)*M*M*M + 2*(Ny + m)*M*M + 2*(Nz + m)*M + 2*(Sz + 1);
}
