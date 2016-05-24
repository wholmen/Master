#include "ccdintermediates.h"

CCDIntermediates::CCDIntermediates()
{

}


CCDIntermediates::CCDIntermediates(basis_set BASIS){
    basis = BASIS;

    // Calculating important variables
    Nholes = basis.Nparticles; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = basis.nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;


    // Weight when adding diagrams to new amplitudes
    weight = 0.5;

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

    // Generate first set of new amplitudes and do second calculation
    UpdateAmplitudes();
    double E1 = CorrelationEnergy();

    // Start the iteration process
    NIterations = 0; tolerance = 1e-6;
    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < MaxIterations){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrelationEnergy();
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

    NIterations = 0; tolerance = 1e-6;
    energies.insert_rows(NIterations,1);
    energies(NIterations) = E1;

    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < 10){

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

                    E += basis.TwoBodyOperator(i,j,a,b) * t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
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

                    // a, b used as states, sent to basis.epsilon and basis.TwoBodyOperator
                    // aa,bb used as iteration elements to pick out state a,b stored in matrices t,t0,I1,I2,I3,I4
                    int a = aa + Nholes; int b = bb + Nholes;
                    double tau = 0;

                    double term = 0;
                    for (int cc=0; cc<Nparticles; cc++){
                        for (int dd=0; dd<Nparticles; dd++){

                            int c = cc + Nholes; int d = dd + Nholes;
                            term += basis.TwoBodyOperator(a,b,c,d) * t0( Index(cc,dd,i,j,Nparticles,Nparticles,Nholes)); // i+j*Nholes, cc+dd*Nparticles);
                        }
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            term += I1(i+j*Nholes, k+l*Nholes) * t0( Index(aa,bb,k,l,Nparticles,Nparticles,Nholes)); //k+l*Nholes, aa+bb*Nparticles);
                        }
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        for (int cc=0; cc<Nparticles; cc++){

                            term += I2(j+k*Nholes, bb+cc*Nparticles) * t0( Index(aa,cc,i,k,Nparticles,Nparticles,Nholes)); //)  i+k*Nholes, aa+cc*Nparticles); // No permutation
                            term -= I2(i+k*Nholes, bb+cc*Nparticles) * t0( Index(aa,cc,j,k,Nparticles,Nparticles,Nholes)); //j+k*Nholes, aa+cc*Nparticles); // Permutation i,j
                            term -= I2(j+k*Nholes, aa+cc*Nparticles) * t0( Index(bb,cc,i,k,Nparticles,Nparticles,Nholes)); //i+k*Nholes, bb+cc*Nparticles); // Permutation a,b
                            term += I2(i+k*Nholes, aa+cc*Nparticles) * t0( Index(bb,cc,j,k,Nparticles,Nparticles,Nholes)); //j+k*Nholes, bb+cc*Nparticles); // Permutation a,b,i,j
                        }
                    }
                    tau += term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        term -= I3(j,k) * t0( Index(aa,bb,i,k,Nparticles,Nparticles,Nholes)); //i+k*Nholes, aa+bb*Nparticles); // No permutation
                        term += I3(i,k) * t0( Index(aa,bb,j,k,Nparticles,Nparticles,Nholes)); //j+k*Nholes, aa+bb*Nparticles); // Permutation i,j
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int cc=0; cc<Nparticles; cc++){
                        term -= I4(bb,cc) * t0( Index(aa,cc,i,j,Nparticles,Nparticles,Nholes)); //) i+j*Nholes, aa+cc*Nparticles); // No permutation
                        term += I4(aa,cc) * t0( Index(bb,cc,i,j,Nparticles,Nparticles,Nholes)); //i+j*Nholes, bb+cc*Nparticles); // Permutation a,b
                    }
                    tau += 0.5*term;

                    tau = basis.TwoBodyOperator(a,b,i,j) + weight*tau; // Weighting the iterative scheme

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes)) = tau / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
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
                            I1(i+j*Nholes, k+l*Nholes) += basis.TwoBodyOperator(k,l,c,d) * t( Index(cc,dd,i,j,Nparticles,Nparticles,Nholes)); //i+j*Nholes, cc+dd*Nparticles);
                        }
                    }
                    I1(i+j*Nholes, k+l*Nholes) = basis.TwoBodyOperator(k,l,i,j) + 0.5*I1(i+j*Nholes, k+l*Nholes);
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
                            I2(j+k*Nholes, bb+cc*Nparticles) += basis.TwoBodyOperator(k,l,c,d) * t( Index(dd,bb,l,j,Nparticles,Nparticles,Nholes)); //l+j*Nholes, dd+bb*Nparticles);
                        }
                    }
                    I2(j+k*Nholes, bb+cc*Nparticles) = basis.TwoBodyOperator(k,b,c,j) + 0.5*I2(j+k*Nholes, bb+cc*Nparticles);
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

                        I3(j,k) += basis.TwoBodyOperator(k,l,c,d) * t( Index(cc,dd,j,l,Nparticles,Nparticles,Nholes)); //j+l*Nholes, cc+dd*Nparticles);
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

                        I4(bb,cc) += basis.TwoBodyOperator(k,l,c,d) * t( Index(bb,dd,k,l,Nparticles,Nparticles,Nholes)); //k+l*Nholes, bb+dd*Nparticles);
                    }
                }
            }
        }
    }
}





double CCDIntermediates::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int CCDIntermediates::Index(int p, int q, int r, int s, int Np, int Nq, int Nr){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}
