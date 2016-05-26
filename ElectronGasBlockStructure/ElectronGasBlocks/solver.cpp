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

    // Weight when adding diagrams to new amplitudes
    weight = 0.5;

    // Setting up two-state configurations

    TwoBodyConfigurations();
    NPARTICLES = Particles.n_rows;
    NHOLES = Holes.n_rows;
    NSTATES = States.n_rows;

    CreateBlocks();

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);
}


double Solver::CCD(bool naive){

    if (naive){
        double E0 = CorrolationEnergy();

        cout << "Energy using blocking: " << E0 << endl;

        UpdateAmplitudes_Naive();
        double E1 = CorrolationEnergy();

        cout << "Energy using blocking: " << E1 << endl;

        NIterations = 0; tolerance = 1e-6;
        while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < 10){

            E0 = E1;
            UpdateAmplitudes_Naive();
            E1 = CorrolationEnergy();
            NIterations ++;

            cout << "Energy using blocking: " << E1 << endl;
        }
        return E1;
    }
    else {
        // Set up the first calculation for all amplitudes equal to zero
        double E0 = CorrolationEnergy(); // Can be hardcoded to 0 to save computation cost

        cout << "Energy using blocking: " << E0 << endl;

        // Generate first set of new amplitudes and do second calculation
        UpdateAmplitudes();
        double E1 = CorrolationEnergy();

        cout << "Energy using blocking: " << E1 << endl;

        // Start the iteration process
        NIterations = 0; tolerance = 1e-6;
        while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < 10){

            E0 = E1;
            UpdateAmplitudes();
            E1 = CorrolationEnergy();
            NIterations ++;

            cout << "Energy using blocking: " << E1 << endl;

        }
        return E1;
    }
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
    CreateBlocksPHPH();
    CreateBlocksHPHP();
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

void Solver::CreateBlocksPHPH(){

    vec identifiers = zeros<vec>(0); Nphph = 0;
    for (int x1 = 0; x1<NX; x1++){
        for (int x2 = 0; x2<NX; x2++){

            if ( Xph(x1,2) == Xph(x2,2) ){ // Identifiers match. They belong to the same block

                bool IDExist = any( identifiers == Xph(x1,2)); // true: Block already added to identifers

                if ( ! IDExist ){ // Block not yet added
                    // Adding block to identifiers
                    identifiers.insert_rows(Nphph,1);
                    identifiers(Nphph) = Xph(x1,2);
                    Nphph ++;
                }
            }
        }
    }
    blocksphph = new Block*[Nphph];

    for (int n=0; n<Nphph; n++){
        blocksphph[n] = new Block(basis, Nholes, Nparticles);
    }

    for (int x1=0; x1<NX; x1++){
        for (int x2=0; x2<NX; x2++){

            if (Xph(x1,2) == Xph(x2,2)) { // Identifiers match. Add states to block

                uvec indices = find( identifiers == Xph(x1,2) );
                int indice = indices(0);

                blocksphph[indice]->AddCrossStates( Xph.row(x1), Xph.row(x2) );
            }
        }
    }
}

void Solver::CreateBlocksHPHP(){

    vec identifiers = zeros<vec>(0); Nhphp = 0;
    for (int x1=0; x1<NX; x1++){
        for (int x2=0; x2<NX; x2++){

            if ( Xhp(x1,2) == Xhp(x2,2) ){

                bool IDExist = any( identifiers == Xhp(x1,2) );

                if ( ! IDExist){

                    identifiers.insert_rows(Nhphp,1);
                    identifiers(Nhphp) = Xhp(x1,2);
                    Nhphp ++;
                }
            }
        }
    }
    blockshphp = new Block*[Nhphp];

    for (int n=0; n<Nhphp; n++){
        blockshphp[n] = new Block(basis, Nholes, Nparticles);
    }
    for (int x1=0; x1<NX; x1++){
        for (int x2=0; x2<NX; x2++){

            if (Xhp(x1,2) == Xhp(x2,2)){

                uvec indices = find ( identifiers == Xhp(x1,2));
                int indice = indices(0);

                blockshphp[indice]->AddCrossStates( Xhp.row(x1), Xhp.row(x2));
            }
        }
    }

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
    }

    // Blocking approach
    for (int n=0; n<Nhphp; n++){
        blockshphp[n]->SetUpMatrices_Lc(t0);
    }
    for (int n=0; n<Nhphp; n++){

        mat Lc = blockshphp[n]->V * blockshphp[n]->T;
        //Lc.print();

    }

}

void Solver::DiagramQa(){

    for (int n=0; n<Npphh; n++){
        blockspphh[n]->SetUpMatrices_Qa(t0);
    }
    // Perform calculations and insert values into t
    for (int n=0; n<Npphh; n++){

        mat Qa = blockspphh[n]->T * blockspphh[n]->V * blockspphh[n]->T2;

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

    // Non-block approach. Testet and working
    mat T1 = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);
    mat V  = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);
    mat T2 = zeros<mat>(Nparticles*Nholes, Nparticles*Nholes);
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    int a = aa + Nholes; int b = bb + Nholes;

                    T1(aa+i*Nparticles, bb+j*Nparticles) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                    V (bb+j*Nparticles, aa+i*Nparticles) = basis.TwoBodyOperator(j,i,b,a);
                    T2(aa+i*Nparticles, bb+j*Nparticles) = t0(Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
                }
            }
        }
    }
    mat Qbt = T1*V*T2;
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qbt( aa+i*Nparticles, bb+j*Nparticles );
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qbt( aa+j*Nparticles, bb+i*Nparticles ); //P(ij)
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) -= weight * 0.5 * Qbt( bb+i*Nparticles, aa+j*Nparticles ); //P(ab)
                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) += weight * 0.5 * Qbt( bb+j*Nparticles, aa+i*Nparticles ); //P(ij)P(ab)
                }
            }
        }
    }

    // Blocking approach
    for (int n=0; n<Nphph; n++){
        blocksphph[n]->SetUpMatrices_Qb(t0);
    }
    for (int n=0; n<Nphph; n++){

        mat Qd = blocksphph[n]->T * blocksphph[n]->V * blocksphph[n]->T2;

        for (double x1=0; x1<blocksphph[n]->X1.n_rows; x1++){
            for (double x2=0; x2<blocksphph[n]->X1.n_rows; x2++){

                int i = blocksphph[n]->X1(x1,1); int j = blocksphph[n]->X1(x2,1);
                int a = blocksphph[n]->X1(x1,0); int b = blocksphph[n]->X1(x2,0); int aa = a-Nholes; int bb = b-Nholes;

                //t( Index(aa,i,bb,j,Nparticles,Nholes,Nparticles)) += weight * 0.5 * Qd(x1,x2) / basis.epsilon(i,j,a,b);
            }
        }


        /*
        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * Qd(A,I) / basis.epsilon(i,j,a,b);
            }
        }*/
    }
}

void Solver::DiagramQc(){

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
    }


    // Blocking approach. Working.
    for (int n=0; n<Npphh; n++){
        blockspphh[n]->SetUpMatrices_Qc(t0);
    }

    for (int n=0; n<Npphh; n++){

        mat Qc_tilde = blockspphh[n]->T * blockspphh[n]->V * blockspphh[n]->T2;

        int Nh = blockspphh[n]->Nh; int Np = blockspphh[n]->Np;

        for (int I=0; I<blockspphh[n]->Nh; I++){
            for (int J=0; J<blockspphh[n]->Nh; J++){
                for (int A=0; A<blockspphh[n]->Np; A++){
                    for (int B=0; B<blockspphh[n]->Np; B++){

                        int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(J,0);
                        int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(B,0);
                        int aa = a-Nholes; int bb = b-Nholes;

                        if ( t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) == 0){
                            //t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) -= weight * 0.5 * Qc_tilde(A + B*Np + I*Np*Np, J);
                            //t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * Qc_tilde(A + B*Np + J*Np*Np, I);
                        }
                    }
                }
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

    // Blocking approach
    for (int n=0; n<Npphh; n++){
        blockspphh[n]->SetUpMatrices_Qd(t0);
    }

    for (int n=0; n<Npphh; n++){

        mat Qd = blockspphh[n]->T * blockspphh[n]->V * blockspphh[n]->T2;

        // Assigning values to t
        for (double I=0; I<blockspphh[n]->Holes.n_rows; I++){
            for (double A=0; A<blockspphh[n]->Particles.n_rows; A++){

                int i = blockspphh[n]->Holes(I,0); int j = blockspphh[n]->Holes(I,1);
                int a = blockspphh[n]->Particles(A,0); int b = blockspphh[n]->Particles(A,1);
                int aa = a-Nholes; int bb = b-Nholes;

                //t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) -= weight * 0.5 * Qd(A,I);
            }
        }
    }
}



void Solver::UpdateAmplitudes_Naive(){
    t0 = t;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    double tau = 0.0;

                    // Adding La
                    for (int cc=0; cc<Nparticles; cc++){
                        for (int dd=0; dd<Nparticles; dd++){
                            int c = cc + Nholes; int d = dd + Nholes;

                            tau += 0.5 * basis.TwoBodyOperator(a,b,c,d) * t0( Index(cc,dd,i,j, Nparticles,Nparticles,Nholes) );
                        }
                    }

                    // Adding Lb
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            tau += 0.5 * basis.TwoBodyOperator(k,l,i,j) * t0( Index(aa,bb,k,l, Nparticles,Nparticles,Nholes) );
                        }
                    }

                    // Adding Lc
                    for (int cc=0; cc<Nparticles; cc++){
                        for (int k=0; k<Nholes; k++){
                            int c = cc + Nholes;

                            tau += basis.TwoBodyOperator(k,b,c,j) * t0( Index(aa,cc,i,k,Nparticles,Nparticles,Nholes)); // No permutation
                            tau -= basis.TwoBodyOperator(k,b,c,i) * t0( Index(aa,cc,j,k,Nparticles,Nparticles,Nholes)); // Permutation of i,j
                            tau -= basis.TwoBodyOperator(k,a,c,j) * t0( Index(bb,cc,i,k,Nparticles,Nparticles,Nholes)); // Permutation of a,b
                            tau += basis.TwoBodyOperator(k,a,c,i) * t0( Index(bb,cc,j,k,Nparticles,Nparticles,Nholes)); // Permutation of a,b,i,j
                        }
                    }


                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){
                            for (int cc=0; cc<Nparticles; cc++){
                                for (int dd=0; dd<Nparticles; dd++){
                                    int c = cc + Nholes; int d = dd + Nholes;

                                    double Qa = 0; double Qb = 0; double Qc = 0; double Qd = 0;

                                    Qa = 0.25*basis.TwoBodyOperator(k,l,c,d) * t0( Index(cc,dd,i,j,Nparticles,Nparticles,Nholes) ) * t0( Index(aa,bb,k,l,Nparticles,Nparticles,Nholes));

                                    Qb += t0( Index(aa,cc,i,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,bb,l,j,Nparticles,Nparticles,Nholes)); // No permutation
                                    Qb -= t0( Index(aa,cc,j,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,bb,l,i,Nparticles,Nparticles,Nholes)); // Permutation of i,j
                                    Qb -= t0( Index(bb,cc,i,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,aa,l,j,Nparticles,Nparticles,Nholes)); // Permutation of a,b
                                    Qb += t0( Index(bb,cc,j,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,aa,l,i,Nparticles,Nparticles,Nholes)); // Permutation of a,b,i,j
                                    Qb *= 0.5 * basis.TwoBodyOperator(k,l,c,d);


                                    Qc -= t0( Index(aa,bb,i,k,Nparticles,Nparticles,Nholes) ) * t0( Index(cc,dd,j,l,Nparticles,Nparticles,Nholes)); // No permutation
                                    Qc += t0( Index(aa,bb,j,k,Nparticles,Nparticles,Nholes) ) * t0( Index(cc,dd,i,l,Nparticles,Nparticles,Nholes)); // Permutation of i,j
                                    Qc *= 0.5 * basis.TwoBodyOperator(k,l,c,d);


                                    Qd -= t0( Index(bb,dd,k,l,Nparticles,Nparticles,Nholes) ) * t0( Index(aa,cc,i,j,Nparticles,Nparticles,Nholes)); // No permutation
                                    Qd += t0( Index(aa,dd,k,l,Nparticles,Nparticles,Nholes) ) * t0( Index(bb,cc,i,j,Nparticles,Nparticles,Nholes)); // Permutation of a,b
                                    Qd *= 0.5 * basis.TwoBodyOperator(k,l,c,d);


                                    tau += Qa + Qb + Qc + Qd;
                                }
                            }
                        }
                    }

                    tau = basis.TwoBodyOperator(a,b,i,j) + weight*tau; //Weighting the iterative scheme

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) = tau / basis.epsilon(i,j,a,b);
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
    States = zeros<mat>(0,3);

    int n=0; // n will count how many two-state combinations we find. Used as indice in the matrix
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            if (i != j){ // Pauli principle demands that the particles must be unequal

                // Setting up direct channels for holes

                // Two-hole momentum and spin
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

    Xhp = zeros<mat>(0,3);
    Xph = zeros<mat>(0,3);

    n=0;
    for (int i=0; i<Nholes; i++){
        for (int aa=0; aa<Nparticles; aa++){
            int a = aa + Nholes;

            // Two-state
            int Nxc = basis.states(i,1) - basis.states(a,1);
            int Nyc = basis.states(i,2) - basis.states(a,2);
            int Nzc = basis.states(i,3) - basis.states(a,3);
            int Szc = basis.states(i,4) - basis.states(a,4);

            Xhp.insert_rows(n,1);
            Xph.insert_rows(n,1);

            Xhp(n,0) = i; Xhp(n,1) = a; Xhp(n,2) = Identifier(Nxc,Nyc,Nzc,Szc);
            Xph(n,0) = a; Xph(n,1) = i; Xph(n,2) = Identifier(-Nxc,-Nyc,-Nzc,-Szc);

            n++;
        }
    }
    NX = Xhp.n_rows;
}

double Solver::Identifier(int Nx, int Ny, int Nz, int Sz){
    return 2*(Nx + m)*M*M*M + 2*(Ny + m)*M*M + 2*(Nz + m)*M + 2*(Sz + 1);
}


























