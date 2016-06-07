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
    NPARTICLES = Particles.n_rows;
    NHOLES = Holes.n_rows;

    CreateBlocks();

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);
}


double Solver::CCD(int MaxIterations){
    // Set up the first calculation for all amplitudes equal to zero
    double E0 = CorrolationEnergy(); // Can be hardcoded to 0 to save computation cost

    cout << "Energy using blocking. E0/N: " << E0 / Nholes << "  E0: " << E0 << endl;

    // Generate first set of new amplitudes and do second calculation
    UpdateAmplitudes();
    double E1 = CorrolationEnergy();

    cout << "Energy using blocking. E0/N: " << E1 / Nholes << "  E0: " << E1 << endl;

    // Start the iteration process
    NIterations = 0; tolerance = 1e-6;
    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < MaxIterations){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrolationEnergy();
        NIterations ++;

        cout << "Energy using blocking. E0/N: " << E1 / Nholes << "  E0: " << E1 << endl;

    }
    return E1;
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

    DiagramL0();

    DiagramLa();
    DiagramLb();
    DiagramLc();

    DiagramQa();
    DiagramQb();
    DiagramQc();
    DiagramQd();

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) /= basis.epsilon(i,j,a,b);
                }
            }
        }
    }
}


void Solver::CreateBlocks(){

    CreateBlocksPPHH();
    CreateBlocksPHHP();
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

                blocksphhp[indice]->AddPHHPStates(Xph.row(x1),Xhp.row(x2));
            }
        }
    }
    for (int n=0; n<Nphhp; n++) blocksphhp[n]->FinishBlock();
}


void Solver::DiagramL0(){


    // Set up matrices for all blockspphh
    for (int i=0; i<Npphh; i++){
        blockspphh[i] ->SetUpMatrices_L0();
    }
    // For L0 diagram, not computation is needed. Only aligning V elements to t.

    //Align elements to t.
    for (int n=0; n<Npphh; n++){

        for (double I=0; I < blockspphh[n]->Holes.n_rows ; I++){
            for (double A=0; A < blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += blockspphh[n]->V(A,I);
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
        mat La = blockspphh[n]->V * blockspphh[n]->T;

        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * La(A,I);
            }
        }
    }
}

void Solver::DiagramLb(){

    for (int i=0; i<Npphh; i++){
        blockspphh[i]->SetUpMatrices_Lb(t0);
    }
    // Do calculation and insert the values into t.
    for (int n=0; n<Npphh; n++){

        mat Lb = blockspphh[n]->T * blockspphh[n]->V;

        // Assigning the values calculated to the right position in t
        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * Lb(A,I);
            }
        }
    }
}

void Solver::DiagramLc(){
    /*
    // Non-block method
    mat T = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);
    mat V = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;

                    T (bb+j*Nparticles, aa+i*Nparticles) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                    V (aa+i*Nparticles, bb+j*Nparticles) = basis.TwoBodyOperator(j,a,b,i);
                }
            }
        }
    }
    mat Lct = V*T;
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * Lct( bb+j*Nparticles, aa+i*Nparticles );
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * Lct( bb+i*Nparticles, aa+j*Nparticles ); //P(ij)
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * Lct( aa+j*Nparticles, bb+i*Nparticles ); //P(ab)
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * Lct( aa+i*Nparticles, bb+j*Nparticles ); //P(ij)P(ab)
                }
            }
        }
    }*/

    // Blocking method
    for (int n=0; n<Nphhp; n++) blocksphhp[n]->SetUpMatrices_Lc(t0);

    for (int n=0; n<Nphhp; n++){

        mat Lc = blocksphhp[n]->V * blocksphhp[n]->T;
        int Nph = blocksphhp[n]->Nph;

        for (int x1=0; x1<Nph; x1++){
            for (int x2=0; x2<Nph; x2++){

                int i=blocksphhp[n]->Xph(x1,1); int j=blocksphhp[n]->Xhp(x2,0);
                int a=blocksphhp[n]->Xph(x1,0); int b=blocksphhp[n]->Xhp(x2,1);
                int aa = a-Nholes; int bb = b-Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * Lc(x1,x2);
                t( Index(aa,bb,j,i,Nparticles,Nparticles,Nholes) ) -= weight * Lc(x1,x2);
                t( Index(bb,aa,i,j,Nparticles,Nparticles,Nholes) ) -= weight * Lc(x1,x2);
                t( Index(bb,aa,j,i,Nparticles,Nparticles,Nholes) ) += weight * Lc(x1,x2);
            }
        }
    }
}

void Solver::DiagramQa(){

    for (int n=0; n<Npphh; n++){
        blockspphh[n]->SetUpMatrices_Qa(t0);
    }
    // Perform calculations and insert values into t
    for (int n=0; n<Npphh; n++){

        mat Qa = blockspphh[n]->T * blockspphh[n]->V * blockspphh[n]->T;

        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.25 * Qa(A,I);
            }
        }
    }
}

void Solver::DiagramQb(){
    /*
    // Non-block approach. Testet and working
    mat T1 = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);
    mat V  = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;

                    T1(aa+i*Nparticles, bb+j*Nparticles) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                    V (bb+j*Nparticles, aa+i*Nparticles) = basis.TwoBodyOperator(j,i,b,a);

                }
            }
        }
    }
    mat Qbt = T1*V*T1;
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    // One can choose wether to implement perturbation perturbing in Qbt or in t(Index)

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qbt( aa+i*Nparticles, bb+j*Nparticles );
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qbt( aa+j*Nparticles, bb+i*Nparticles ); //P(ij)
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qbt( bb+i*Nparticles, aa+j*Nparticles ); //P(ab)
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qbt( bb+j*Nparticles, aa+i*Nparticles ); //P(ij)P(ab)

                    //t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qbt( aa+i*Nparticles, bb+j*Nparticles );
                    //t( Index(aa,bb,j,i,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qbt( aa+i*Nparticles, bb+j*Nparticles ); //P(ij)
                    //t( Index(bb,aa,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qbt( aa+i*Nparticles, bb+j*Nparticles ); //P(ab)
                    //t( Index(bb,aa,j,i,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qbt( aa+i*Nparticles, bb+j*Nparticles ); //P(ij)P(ab)

                }
            }
        }
    }*/

    // Blocking approach
    for (int n=0; n<Nphhp; n++) blocksphhp[n]->SetUpMatrices_Qb(t0);

    for (int n=0; n<Nphhp; n++){

        mat Qb = blocksphhp[n]->T * blocksphhp[n]->V * blocksphhp[n]->T;
        int Nph = blocksphhp[n]->Nph;

        for (int x1=0; x1<Nph; x1++){
            for (int x2=0; x2<Nph; x2++){

                int i=blocksphhp[n]->Xph(x1,1); int j=blocksphhp[n]->Xhp(x2,0);
                int a=blocksphhp[n]->Xph(x1,0); int b=blocksphhp[n]->Xhp(x2,1);
                int aa = a-Nholes; int bb = b-Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * Qb(x1,x2);
                t( Index(aa,bb,j,i,Nparticles,Nparticles,Nholes) ) -= weight * 0.5 * Qb(x1,x2);
                t( Index(bb,aa,i,j,Nparticles,Nparticles,Nholes) ) -= weight * 0.5 * Qb(x1,x2);
                t( Index(bb,aa,j,i,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * Qb(x1,x2);
            }
        }
    }
}

void Solver::DiagramQc(){
    /*
    // Non-blocking approach. Tested and working
    mat T1 = zeros<mat>(Nparticles*Nparticles*Nholes, Nholes);
    mat V  = zeros<mat>(Nholes, Nparticles*Nparticles*Nholes);
    mat T2 = zeros<mat>(Nparticles*Nparticles*Nholes, Nholes);
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;

                    T1(aa + bb*Nparticles + i*Nparticles2, j) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                    V (i, aa + bb*Nparticles + j*Nparticles2) = basis.TwoBodyOperator(i,j,a,b);
                    T2(aa + bb*Nparticles + j*Nparticles2, i) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                }
            }
        }
    }
    mat Qct = T1*V*T2;
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qct( aa+bb*Nparticles+i*Nparticles2, j );
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qct( aa+bb*Nparticles+j*Nparticles2, i );
                }
            }
        }
    }*/

    // Block Approach

    vec identifiers = zeros<vec>(0); Nkhpp = 0;
    for (int k1=0; k1<NK; k1++){
        for (int k2=0; k2<NK3; k2++){

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

    for (int k1=0; k1<NK; k1++){
        for (int k2=0; k2<NK3; k2++){

            if (Kh(k1,1) == Khpp(k2,3) ) { // Block found

                uvec indices = find( identifiers == Kh(k1,1));
                int indice = indices(0);

                blockskhpp[indice]->AddTripleStates( Kh.row(k1), Khpp.row(k2) );
            }
        }
    }
    for (int n=0; n<Nkhpp; n++) blockskhpp[n]->FinishBlock();

    for (int n=0; n<Nkhpp; n++) blockskhpp[n]->SetUpMatrices_Qc(t0);

    for (int n=0; n<Nkhpp; n++){

        mat Qc = blockskhpp[n]->T * blockskhpp[n]->V * blockskhpp[n]->T;

        int Nkh = blockskhpp[n]->Nkh; int Nkhpp = blockskhpp[n]->Nkhpp;

        for (int k1=0; k1<Nkh; k1++){
            for (int k2=0; k2<Nkhpp; k2++){

                int i=blockskhpp[n]->Kh(k1,0);   int j=blockskhpp[n]->Khpp(k2,0);
                int a=blockskhpp[n]->Khpp(k2,1); int b=blockskhpp[n]->Khpp(k2,2);
                int aa=a-Nholes;  int bb=b-Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qc( k2,k1 );
                t( Index(aa,bb,j,i,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qc( k2,k1 );

            }
        }

    }
}

void Solver::DiagramQd(){

    // Non-blocking approach.
    mat T1 = zeros<mat>(Nparticles, Nholes2*Nparticles);
    mat T2 = zeros<mat>(Nparticles, Nholes2*Nparticles);
    mat  V = zeros<mat>(Nholes2*Nparticles, Nparticles);
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;

                    T1(aa, i+j*Nholes+bb*Nholes2) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                    T2(bb, i+j*Nholes+aa*Nholes2) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                    V (i+j*Nholes+bb*Nholes2, aa) = basis.TwoBodyOperator(i,j,a,b);
                }
            }
        }
    }
    mat Qdt = T1*V*T2;
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qdt( bb, i+j*Nholes+aa*Nholes2);
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qdt( aa, i+j*Nholes+bb*Nholes2);
                }
            }
        }
    }
}




double Solver::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int Solver::Index(int p, int q, int r, int s, int Np, int Nq, int Nr){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
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

    Kh = zeros<mat>(0,2);
    Khpp = zeros<mat>(0,4);

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
    NK3 = Khpp.n_rows; NK = Kh.n_rows;
}

double Solver::Identifier(int Nx, int Ny, int Nz, int Sz){
    return 2*(Nx + m)*M*M*M + 2*(Ny + m)*M*M + 2*(Nz + m)*M + 2*(Sz + 1);
}


























