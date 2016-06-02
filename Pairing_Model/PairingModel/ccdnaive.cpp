#include "ccdnaive.h"

CCDNaive::CCDNaive(){}

CCDNaive::CCDNaive(basis_set BASIS){
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
}

vec CCDNaive::CCD_ReturnAllIterations(){
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

double CCDNaive::CCD(int MaxIterations){
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


double CCDNaive::CorrelationEnergy(){
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


void CCDNaive::UpdateAmplitudes(){
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

                    // Adding quadratic terms
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

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) = tau / basis.epsilonijab(i,j,a,b);
                }
            }
        }
    }
}

double CCDNaive::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int CCDNaive::Index(int p, int q, int r, int s, int Np, int Nq, int Nr){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}
