#include "ccdintermediates.h"

CCDIntermediates::CCDIntermediates(){}

CCDIntermediates::CCDIntermediates(PairingBasis BASIS){
    pabasis = BASIS; BasisNumber = 1;

    // Calculating important variables
    Nholes = pabasis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = pabasis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

    // Weight when adding diagrams to new amplitudes
    weight = 1.0;

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);

    // Setting up matrices for intermediate calculation
    I1 = zeros<mat>(Nholes*Nholes, Nholes*Nholes);
    I2 = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);
    I3 = zeros<mat>(Nholes, Nholes);
    I4 = zeros<mat>(Nparticles, Nparticles);
}

CCDIntermediates::CCDIntermediates(ElectronBasis BASIS){
    elbasis = BASIS; BasisNumber = 2;

    // Calculating important variables
    Nholes = elbasis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = elbasis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

    // Weight when adding diagrams to new amplitudes
    weight = 1.0;

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);

    // Setting up matrices for intermediate calculation
    I1 = zeros<mat>(Nholes*Nholes, Nholes*Nholes);
    I2 = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);
    I3 = zeros<mat>(Nholes, Nholes);
    I4 = zeros<mat>(Nparticles, Nparticles);
}

CCDIntermediates::CCDIntermediates(NuclearBasis BASIS){
    ncbasis = BASIS; BasisNumber = 3;

    // Calculating important variables
    Nholes = ncbasis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = ncbasis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

    // Weight when adding diagrams to new amplitudes
    weight = 1.0;

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);

    // Setting up matrices for intermediate calculation
    I1 = zeros<mat>(Nholes*Nholes, Nholes*Nholes);
    I2 = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);
    I3 = zeros<mat>(Nholes, Nholes);
    I4 = zeros<mat>(Nparticles, Nparticles);
}


double CCDIntermediates::CCD(int MaxIterations){
    // Set up the first calculation for all amplitudes equal to zero
    double E0 = CorrelationEnergy(); // Can be hardcoded to 0 to save computation cost
    cout << "Energy using intermediates. E0/N: " << E0 / Nholes << "  E0: " << E0 << endl;

    // Generate first set of new amplitudes and do second calculation
    UpdateAmplitudes();
    double E1 = CorrelationEnergy();
    cout << "Energy using intermediates. E0/N: " << E1 / Nholes << "  E0: " << E1 << endl;

    // Start the iteration process
    NIterations = 0; tolerance = 1e-12;
    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < MaxIterations){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrelationEnergy();
        cout << "Energy using intermediates. E0/N: " << E1 / Nholes << "  E0: " << E1 << endl;
        NIterations ++;
    }
    return E1;
}

vec CCDIntermediates::CCD_ReturnAllIterations(){
    // Variation of CCD() made to store correlation energy for every iteration
    vec energies = zeros<vec>(0);

    double E0 = CorrelationEnergy();
    UpdateAmplitudes();
    double E1 = CorrelationEnergy();

    NIterations = 0; tolerance = 1e-12;
    energies.insert_rows(NIterations,1);
    energies(NIterations) = E1; NIterations++;

    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < 20){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrelationEnergy();
        energies.insert_rows(NIterations,1);
        energies(NIterations) = E1;

        NIterations ++;
    }
    return energies;
}

double CCDIntermediates::CorrelationEnergy(){
    double E = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += v(i,j,a,b) * t( Index(aa,bb,i,j));
                }
            }
        }
    }
    return E * 0.25;
}

void CCDIntermediates::UpdateAmplitudes(){

    t0 = t; UpdateI1(); UpdateI2(); UpdateI3(); UpdateI4();

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){

                    // a, b used as states, sent to basis.epsilon and v
                    // aa,bb used as iteration elements to pick out state a,b stored in matrices t,t0,I1,I2,I3,I4
                    int a = aa + Nholes; int b = bb + Nholes;
                    double tau = 0;

                    double term = 0;
                    for (int cc=0; cc<Nparticles; cc++){
                        for (int dd=0; dd<Nparticles; dd++){

                            int c = cc + Nholes; int d = dd + Nholes;
                            term += v(a,b,c,d) * t0( Index(cc,dd,i,j)); // i+j*Nholes, cc+dd*Nparticles);
                        }
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            term += I1(i+j*Nholes, k+l*Nholes) * t0( Index(aa,bb,k,l)); //k+l*Nholes, aa+bb*Nparticles);
                        }
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        for (int cc=0; cc<Nparticles; cc++){

                            term += I2(j+k*Nholes, bb+cc*Nparticles) * t0( Index(aa,cc,i,k)); //)  i+k*Nholes, aa+cc*Nparticles); // No permutation
                            term -= I2(i+k*Nholes, bb+cc*Nparticles) * t0( Index(aa,cc,j,k)); //j+k*Nholes, aa+cc*Nparticles); // Permutation i,j
                            term -= I2(j+k*Nholes, aa+cc*Nparticles) * t0( Index(bb,cc,i,k)); //i+k*Nholes, bb+cc*Nparticles); // Permutation a,b
                            term += I2(i+k*Nholes, aa+cc*Nparticles) * t0( Index(bb,cc,j,k)); //j+k*Nholes, bb+cc*Nparticles); // Permutation a,b,i,j
                        }
                    }
                    tau += term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        term -= I3(j,k) * t0( Index(aa,bb,i,k)); //i+k*Nholes, aa+bb*Nparticles); // No permutation
                        term += I3(i,k) * t0( Index(aa,bb,j,k)); //j+k*Nholes, aa+bb*Nparticles); // Permutation i,j
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int cc=0; cc<Nparticles; cc++){
                        term -= I4(bb,cc) * t0( Index(aa,cc,i,j)); //) i+j*Nholes, aa+cc*Nparticles); // No permutation
                        term += I4(aa,cc) * t0( Index(bb,cc,i,j)); //i+j*Nholes, bb+cc*Nparticles); // Permutation a,b
                    }
                    tau += 0.5*term;

                    tau += v(a,b,i,j); // Weighting the iterative scheme

                    t( Index(aa,bb,i,j)) = tau / epsilon(i,j,a,b);
                }
            }
        }
    }
    // Dividing by factor epsilon(i,j,a,b)
    //t = t / epsilon();

    // Adding weight factor
    if (weight != 0) t = weight*t + (1-weight)*t0;
}

void CCDIntermediates::UpdateI1(){
    // I1 is the Intermediate matrix 1. It contains pre calculated values for all variations of
    // i, j, k, l. The values are located at I1(i + j*Nholes, k + l*Nholes)
    // I1 has the size (Nholes^2, Nholes^2)
    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int k=0; k<Nholes; k++){
                for (int l=0; l<Nholes; l++){

                    I1(i+j*Nholes, k+l*Nholes) = 0;

                    for (int cc=0; cc<Nparticles; cc++){
                        for (int dd=0; dd<Nparticles; dd++){

                            int c = cc + Nholes; int d = dd + Nholes;
                            I1(i+j*Nholes, k+l*Nholes) += v(k,l,c,d) * t( Index(cc,dd,i,j)); //i+j*Nholes, cc+dd*Nparticles);
                        }
                    }
                    I1(i+j*Nholes, k+l*Nholes) = v(k,l,i,j) + 0.5*I1(i+j*Nholes, k+l*Nholes);
                }
            }
        }
    }
}

void CCDIntermediates::UpdateI2(){
    // I2 is the Intermediate matrix 2. It contains pre-calculated values for all variations of
    // j, k, b, c. The values are located at I2(j + k*Nholes, b + c*Nparticles)
    // I2 has the size (Nholes^2, Nparticles^2)
    for (int j=0; j<Nholes; j++){
        for (int k=0; k<Nholes; k++){
            for (int bb=0; bb<Nparticles; bb++){
                for (int cc=0; cc<Nparticles; cc++){
                    int b = bb + Nholes; int c = cc + Nholes; // Converting iteration element into basis state number

                    I2(j+k*Nholes, bb+cc*Nparticles) = 0;

                    for (int l=0; l<Nholes; l++){
                        for (int dd=0; dd<Nparticles; dd++){

                            int d = dd + Nholes;
                            I2(j+k*Nholes, bb+cc*Nparticles) += v(k,l,c,d) * t( Index(dd,bb,l,j)); //l+j*Nholes, dd+bb*Nparticles);
                        }
                    }
                    I2(j+k*Nholes, bb+cc*Nparticles) = v(k,b,c,j) + 0.5*I2(j+k*Nholes, bb+cc*Nparticles);
                }
            }
        }
    }
}

void CCDIntermediates::UpdateI3(){
    // I3 is the Intermediate matrix 2. It contains pre-calculated values for all variations of
    // j, k. The values are located at I2(j,k)
    // I3 has the size (Nholes, Nholes)

    for (int j=0; j<Nholes; j++){
        for (int k=0; k<Nholes; k++){
            I3(j,k) = 0;

            for (int l=0; l<Nholes; l++){
                for (int cc=0; cc<Nparticles; cc++){
                    for (int dd=0; dd<Nparticles; dd++){

                        int c = cc + Nholes; int d = dd + Nholes;

                        I3(j,k) += v(k,l,c,d) * t( Index(cc,dd,j,l)); //j+l*Nholes, cc+dd*Nparticles);
                    }
                }
            }
        }
    }
}

void CCDIntermediates::UpdateI4(){
    // I4 is the Intermediate matrix 2. It contains pre-calculated values for all variations of
    // b, c. The values are located at I4(bb,cc)
    // I4 has the size (Nparticles, Nparticles)

    for (int bb=0; bb<Nparticles; bb++){
        for (int cc=0; cc<Nparticles; cc++){

            I4(bb,cc) = 0;

            for (int k=0; k<Nholes; k++){
                for (int l=0; l<Nholes; l++){
                    for (int dd=0; dd<Nparticles; dd++){

                        int c = cc+Nholes; double d = dd+Nholes;

                        I4(bb,cc) += v(k,l,c,d) * t( Index(bb,dd,k,l)); //k+l*Nholes, bb+dd*Nparticles);
                    }
                }
            }
        }
    }
}

double CCDIntermediates::epsilon(int i, int j, int a, int b){
    if (BasisNumber == 1) return pabasis.epsilon(i,j,a,b);
    else if (BasisNumber == 2) return elbasis.epsilon(i,j,a,b);
    else if (BasisNumber == 3) return ncbasis.epsilon(i,j,a,b);
    else {cout << "basis is not defined properly in ccd naive. Epsilon not computed properly" << endl; return 0;}
}

vec CCDIntermediates::epsilon(){

    if (BasisNumber == 1) return pabasis.EpsilonMatrix;
    else if (BasisNumber == 2) return elbasis.EpsilonMatrix;
    else if (BasisNumber == 3) return ncbasis.EpsilonMatrix;
    else {cout << "basis is not defined properly in ccd naive. Epsilon not computed properly" << endl; return 0;}
}

double CCDIntermediates::v(int p, int q, int r, int s){

    if (BasisNumber == 1) return pabasis.TwoBodyOperator(p,q,r,s);
    else if (BasisNumber == 2) return elbasis.TwoBodyOperator(p,q,r,s);
    else if (BasisNumber == 3) return ncbasis.TwoBodyOperator(p,q,r,s);
    else {cout << "basis is not defined properly in ccd naive. TwoBodyOperator not computed properly" << endl; return 0;}
}

double CCDIntermediates::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int CCDIntermediates::Index(int p, int q, int r, int s){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Nparticles + r*Nparticles2 + s*Nparticles2*Nholes;
}
