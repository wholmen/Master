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
    for (int i=0; i<NidDE; i++){
        blocksDE[i] ->SetUpMatrices_Energy(t);
    }

    // Do the matrix-matrix multiplications for all blocks
    double E = 0;
    for (int i=0; i<NidDE; i++){

        mat Energy = blocksDE[i]->V * blocksDE[i]->T;

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
    DiagramQa();
    DiagramQc();
}



void Solver::CreateBlocks(){
    CreateBlocksDE();
    CreateBlocksL0();
    CreateBlocksLa();
    CreateBlocksLb();
    CreateBlocksQa();
    CreateBlocksQc();
}

void Solver::CreateBlocksDE(){

    vec identifiers = zeros<vec>(0); NidDE = 0;
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) { // They have matching identifier and belong to the same block.

                bool IDExist = any( identifiers == Holes(I,2) ); // Checking wether the identifier has already been counted

                if ( ! IDExist ){ // We have found a new identifier. Storing it in the list of identifiers

                    identifiers.insert_rows(NidDE, 1);
                    identifiers(NidDE) = Holes(I,2);
                    NidDE ++;
                }
            }
        }
    }
    blocksDE = new Block*[NidDE];

    // Set up blocks
    for (int n=0; n<NidDE; n++){
        blocksDE[n] = new Block(basis,Nholes,Nparticles);
    }

    // Add states to blocks
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) {

                uvec indices = find( identifiers == Holes(I,2) );
                int indice = indices(0);

                blocksDE[indice] ->AddStates( Holes.row(I), Particles.row(A) );
            }
        }
    }
}

void Solver::CreateBlocksL0(){

    // Counting and storing all unique identifiers that exist for both Holes and Particles
    vec identifiers = zeros<vec>(0); NidL0 = 0;
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) { // They have matching identifier and belong to the same block.

                bool IDExist = any( identifiers == Holes(I,2) ); // Checking wether the identifier has already been counted

                if ( ! IDExist ){ // We have found a new identifier. Storing it in the list of identifiers

                    identifiers.insert_rows(NidL0, 1);
                    identifiers(NidL0) = Holes(I,2);
                    NidL0 ++;
                }
            }
        }
    }
    blocksL0 = new Block*[NidL0];

    // Initiate all the blocks needed
    for (int i=0; i<NidL0; i++){
        blocksL0[i] = new Block(basis, Nholes, Nparticles);
    }
    // Add all two-state configurations to the corresponding blocksL0
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) {

                uvec indices = find( identifiers == Holes(I,2) );
                int indice = indices(0);

                blocksL0[indice] ->AddStates( Holes.row(I), Particles.row(A) );
            }
        }
    }
}

void Solver::CreateBlocksLa(){

    vec identifiers = zeros<vec>(0); NidLa = 0;
    for (int A=0; A<NPARTICLES; A++){
        for (int C=0; C<NPARTICLES; C++){
            for (int I=0; I<NHOLES; I++){

                if ( Particles(A,2) == Particles(C,2) && Particles(A,2) == Holes(I,2) ){ // Matching identifiers. Must be set up as block

                    bool IDExist = any( identifiers == Holes(I,2) ); // Checking wether the identifier has already been counted

                    if ( ! IDExist ){ // We have found a new identifier. Storing it in the list of identifiers

                        identifiers.insert_rows(NidLa, 1);
                        identifiers(NidLa) = Holes(I,2);
                        NidLa ++;
                    }
                }
            }
        }
    }
    blocksLa = new Block*[NidLa];

    // Initiate all the blocksLa needed
    for (int i=0; i<NidLa; i++){
        blocksLa[i] = new Block(basis, Nholes, Nparticles);
    }
    // Adding states to the blocks
    for (int A=0; A<NPARTICLES; A++){
        for (int I=0; I<NHOLES; I++){
            for (int C=0; C<NPARTICLES; C++){

                if ( Particles(A,2) == Particles(C,2) && Particles(A,2) == Holes(I,2) ){

                    uvec indices = find( identifiers == Particles(A,2) );
                    int indice = indices(0);

                    blocksLa[indice]->AddStates( Holes.row(I), Particles.row(A) );
                }
            }
        }
    }
}

void Solver::CreateBlocksLb(){

    vec identifiers = zeros<vec>(0); NidLb = 0;
    for (int A=0; A<NPARTICLES; A++){
        for (int I=0; I<NHOLES; I++){
            for (int K=0; K<NHOLES; K++){

                if ( Holes(I,2) == Particles(A,2) && Holes(I,2) == Holes(K,2) ){ // They belong to the same block

                    bool IDExist = any( identifiers == Holes(I,2) ); // Has this identifier already been stored?

                    if ( ! IDExist ){ // Identifier is not yet stored. Must be added
                        identifiers.insert_rows(NidLb, 1);
                        identifiers(NidLb) = Holes(I,2); // Saving the identifier for this block
                        NidLb ++;
                    }
                }
            }
        }
    }
    blocksLb = new Block*[NidLb]; // List of all the blocks for diagram Lb

    // Initiate the blocks
    for (int i=0; i<NidLb; i++){
        blocksLb[i] = new Block(basis,Nholes,Nparticles);
    }
    // Adding states to the blocks
    for (int I=0; I<NHOLES; I++){
        for (int K=0; K<NHOLES; K++){
            for (int A=0; A<NPARTICLES; A++){

                if ( Holes(I,2) == Particles(A,2) && Holes(I,2) == Holes(K,2) ){ // They belong to the same block

                    uvec indices = find( identifiers == Holes(I,2) );
                    int indice = indices(0);

                    blocksLb[indice]->AddStates( Holes.row(I), Particles.row(A) );
                }
            }
        }
    }
}

void Solver::CreateBlocksQa(){
    /*
    vec identifiers = zeros<vec>(0); NidQa = 0;
    for (int A=0; A<NPARTICLES; A++){
        for (int C=0; C<NPARTICLES; C++){
            for (int I=0; I<NHOLES; I++){
                for (int K=0; K<NHOLES; K++){

                    if ( Holes(I,2) == Particles(A,2) && Holes(I,2) == Particles(C,2) && Holes(I,2) == Holes(K,2) ){ // I,A,C,K belong to the same block

                        bool IDExist = any( identifiers == Holes(I,2) ); // Is this block already stored?

                        if ( ! IDExist ){
                            identifiers.insert_rows(NidQa, 1);
                            identifiers(NidQa) = Holes(I,2);
                            NidQa ++;
                        }
                    }
                }
            }
        }
    }*/
    vec identifiers = zeros<vec>(0); NidQa = 0;
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) { // They have matching identifier and belong to the same block.

                bool IDExist = any( identifiers == Holes(I,2) ); // Checking wether the identifier has already been counted

                if ( ! IDExist ){ // We have found a new identifier. Storing it in the list of identifiers

                    identifiers.insert_rows(NidQa, 1);
                    identifiers(NidQa) = Holes(I,2);
                    NidQa ++;
                }
            }
        }
    }
    blocksQa = new Block*[NidQa];

    // Initiate the blocks
    for (int n=0; n<NidQa; n++){
        blocksQa[n] = new Block(basis,Nholes,Nparticles);
    }
    /*
    // Add states to blocks
    for (int A=0; A<NPARTICLES; A++){
        for (int C=0; C<NPARTICLES; C++){
            for (int I=0; I<NHOLES; I++){
                for (int K=0; K<NHOLES; K++){

                    if ( Holes(I,2) == Particles(A,2) && Holes(I,2) == Particles(C,2) && Holes(I,2) == Holes(K,2)){

                        uvec indices = find( identifiers == Holes(I,2) );
                        int indice = indices(0);

                        blocksQa[indice]->AddStates( Holes.row(I), Particles.row(A) );
                    }
                }
            }
        }
    }*/
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if ( Holes(I,2) == Particles(A,2) ) {

                uvec indices = find( identifiers == Holes(I,2) );
                int indice = indices(0);

                blocksQa[indice]->AddStates( Holes.row(I), Particles.row(A) );
            }
        }
    }
}

void Solver::CreateBlocksQc(){

    vec identifiers = zeros<vec>(0); NidQc = 0;
    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if (Holes(I,2) == Particles(A,2) ) { // Block found

                bool IDExist = any( identifiers == Holes(I,2) );

                if ( ! IDExist){
                    identifiers.insert_rows(NidQc,1);
                    identifiers(NidQc) = Holes(I,2);
                    NidQc ++;
                }
            }
        }
    }
    blocksQc = new Block*[NidQc];

    for (int n=0; n<NidQc; n++){
        blocksQc[n] = new Block(basis, Nholes, Nparticles);
    }

    for (int I=0; I<NHOLES; I++){
        for (int A=0; A<NPARTICLES; A++){

            if (Holes(I,2) == Particles(A,2)){

                uvec indices = find( identifiers == Holes(I,2) );
                int indice = indices(0);

                blocksQc[indice]->AddStates( Holes.row(I), Particles.row(I) );
            }
        }
    }
}


void Solver::DiagramL0(){


    // Set up matrices for all blocksL0
    for (int i=0; i<NidL0; i++){
        blocksL0[i] ->SetUpMatrices_L0();
    }
    // For L0 diagram, not computation is needed. Only aligning V elements to t.

    //Align elements to t.
    for (int n=0; n<NidL0; n++){

        for (double I=0; I < blocksL0[n]->Holes.n_rows ; I++){
            for (double A=0; A < blocksL0[n]->Particles.n_rows; A++){

                int i = blocksL0[n]->Holes(I,0); int j = blocksL0[n]->Holes(I,1);
                int a = blocksL0[n]->Particles(A,0); int b = blocksL0[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += blocksL0[n]->V(A,I) / basis.epsilon(i,j,a,b);
            }
        }
    }
}

void Solver::DiagramLa(){

    // Tested. Working.

    for (int i=0; i<NidLa; i++){
        blocksLa[i]->SetUpMatrices_La(t0);
    }

    for (int n=0; n<NidLa; n++){
        mat La = blocksLa[n]->V * blocksLa[n]->T;

        for (double I=0; I<blocksLa[n]->Holes.n_rows; I++){
            for (double A=0; A<blocksLa[n]->Particles.n_rows; A++){

                int i = blocksLa[n]->Holes(I,0); int j = blocksLa[n]->Holes(I,1);
                int a = blocksLa[n]->Particles(A,0); int b = blocksLa[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * La(A,I) / basis.epsilon(i,j,a,b);
            }
        }
    }
}

void Solver::DiagramLb(){

    for (int i=0; i<NidLb; i++){
        blocksLb[i]->SetUpMatrices_Lb(t0);
    }
    // Do calculation and insert the values into t.
    for (int n=0; n<NidLb; n++){

        mat LbT = blocksLb[n]->V.t() * blocksLb[n]->T.t(); // Calculating the product of the transposed matrices V.t and T.t
        mat Lb = LbT.t(); // Transposing back to get the diagram Lb

        // Assigning the values calculated to the right position in t
        for (double I=0; I<blocksLb[n]->Holes.n_rows; I++){
            for (double A=0; A<blocksLb[n]->Particles.n_rows; A++){

                int i = blocksLb[n]->Holes(I,0); int j = blocksLb[n]->Holes(I,1);
                int a = blocksLb[n]->Particles(A,0); int b = blocksLb[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.5 * Lb(A,I) / basis.epsilon(i,j,a,b);
            }
        }
    }
}

void Solver::DiagramQa(){

    for (int n=0; n<NidQa; n++){
        blocksQa[n]->SetUpMatrices_Qa(t0);
    }
    // Perform calculations and insert values into t
    for (int n=0; n<NidQa; n++){

        mat Qa = blocksQa[n]->T * blocksQa[n]->V * blocksQa[n]->T2;

        for (double I=0; I<blocksQa[n]->Holes.n_rows; I++){
            for (double A=0; A<blocksQa[n]->Particles.n_rows; A++){

                int i = blocksQa[n]->Holes(I,0); int j = blocksQa[n]->Holes(I,1);
                int a = blocksQa[n]->Particles(A,0); int b = blocksQa[n]->Particles(A,1);
                int aa = a - Nholes; int bb = b - Nholes;

                t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) += weight * 0.25 * Qa(A,I) / basis.epsilon(i,j,a,b);
            }
        }
    }
}

void Solver::DiagramQc(){

}

void Solver::UpdateAmplitudes_Naive(){
    t0 = t;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    double tau = 0.0;
                    /*
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
                    } */
/*
                    // Adding Lc
                    for (int cc=0; cc<Nparticles; cc++){
                        for (int k=0; k<Nholes; k++){
                            int c = cc + Nholes;

                            tau += basis.TwoBodyOperator(k,b,c,j)*t0(i+k*Nholes, aa+cc*Nparticles); // No permutation
                            tau -= basis.TwoBodyOperator(k,b,c,i)*t0(j+k*Nholes, aa+cc*Nparticles); // Permutation of i,j
                            tau -= basis.TwoBodyOperator(k,a,c,j)*t0(i+k*Nholes, bb+cc*Nparticles); // Permutation of a,b
                            tau += basis.TwoBodyOperator(k,a,c,i)*t0(j+k*Nholes, bb+cc*Nparticles); // Permutation of a,b,i,j

                        }
                    }
*/
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){
                            for (int cc=0; cc<Nparticles; cc++){
                                for (int dd=0; dd<Nparticles; dd++){
                                    int c = cc + Nholes; int d = dd + Nholes;

                                    //double Qa = 0; double Qb = 0; double Qc = 0; double Qd = 0;
                                    double Qc = 0;
                                    //Qa = 0.25*basis.TwoBodyOperator(k,l,c,d) * t0( Index(cc,dd,i,j,Nparticles,Nparticles,Nholes) ) * t0( Index(aa,bb,k,l,Nparticles,Nparticles,Nholes));
                                    /*
                                    Qb += t0(i+k*Nholes, aa+cc*Nparticles) * t0(l+j*Nholes, dd+bb*Nparticles); // No permutation
                                    Qb -= t0(j+k*Nholes, aa+cc*Nparticles) * t0(l+i*Nholes, dd+bb*Nparticles); // Permutation of i,j
                                    Qb -= t0(i+k*Nholes, bb+cc*Nparticles) * t0(l+j*Nholes, dd+aa*Nparticles); // Permutation of a,b
                                    Qb += t0(j+k*Nholes, bb+cc*Nparticles) * t0(l+i*Nholes, dd+aa*Nparticles); // Permutation of a,b,i,j
                                    Qb *= 0.5 * basis.TwoBodyOperator(k,l,c,d);
                                    */
                                    Qc -= t0( Index(aa,bb,i,k,Nparticles,Nparticles,Nholes) ) * t0( Index(cc,dd,j,l,Nparticles,Nparticles,Nholes)); // No permutation
                                    //Qc += t0(j+k*Nholes, aa+bb*Nparticles) * t0(i+l*Nholes, cc+dd*Nparticles); // Permutation of i,j
                                    Qc *= 0.5 * basis.TwoBodyOperator(k,l,c,d);
                                    /*
                                    Qd -= t0(k+l*Nholes, bb+dd*Nparticles) * t0(i+j*Nholes, aa+cc*Nparticles); // No permutation
                                    Qd += t0(k+l*Nholes, aa+dd*Nparticles) * t0(i+j*Nholes, bb+cc*Nparticles); // Permutation of a,b
                                    Qd *= 0.5 * basis.TwoBodyOperator(k,l,c,d);
                                    */
                                    tau += Qc; //Qa + Qb + Qc + Qd;
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


























