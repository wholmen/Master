#include "ccdnaive.h"

CCDNaive::CCDNaive(){}

CCDNaive::CCDNaive(PairingBasis BASIS) {
    pabasis = BASIS; BasisNumber = 1;

    // Calculating important variables
    Nholes = pabasis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = pabasis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;


    // Weight when adding diagrams to new amplitudes
    weight = 1.0; tolerance = 1e-6;

    // Setting up matrices
    t0 = zeros<vec>(Nparticles2*Nholes2);
    t  = zeros<vec>(Nparticles2*Nholes2);
}
CCDNaive::CCDNaive(ElectronBasis BASIS){
    elbasis = BASIS; BasisNumber = 2;

    // Calculating important variables
    Nholes = elbasis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = elbasis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

    // Weight when adding diagrams to new amplitudes
    weight = 1.0; tolerance = 1e-6;

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

double CCDNaive::CCD(int MaxIterations){
    // Set up the first calculation for all amplitudes equal to zero
    double E0 = CorrelationEnergy(); // Can be hardcoded to 0 to save computation cost
    cout << "Energy using naive. E0/N: " << E0 / Nholes << "  E0: " << E0 << endl;

    // Generate first set of new amplitudes and do second calculation
    UpdateAmplitudes();
    double E1 = CorrelationEnergy(); //cout << E1 << endl;
    cout << "Energy using naive. E0/N: " << E1 / Nholes << "  E0: " << E1 << endl;

    // Start the iteration process
    NIterations = 0;
    while ( AbsoluteDifference(E1,E0) > tolerance && NIterations < MaxIterations){

        E0 = E1;
        UpdateAmplitudes();
        E1 = CorrelationEnergy();
        cout << "Energy using naive. E0/N: " << E1 / Nholes << "  E0: " << E1 << endl;
        NIterations ++;
    }
    return E1;
}


double CCDNaive::CorrelationEnergy(){
    double E = 0.0;

    if (BasisNumber == 1) PairingBasis basis = pabasis;
    else ElectronBasis basis = elbasis;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += v(i,j,a,b) * t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes));
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

                            tau += 0.5 * v(a,b,c,d) * t0( Index(cc,dd,i,j, Nparticles,Nparticles,Nholes) );
                        }
                    }

                    // Adding Lb
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            tau += 0.5 * v(k,l,i,j) * t0( Index(aa,bb,k,l, Nparticles,Nparticles,Nholes) );
                        }
                    }

                    // Adding Lc
                    for (int cc=0; cc<Nparticles; cc++){
                        for (int k=0; k<Nholes; k++){
                            int c = cc + Nholes;

                            tau += v(k,b,c,j) * t0( Index(aa,cc,i,k,Nparticles,Nparticles,Nholes)); // No permutation
                            tau -= v(k,b,c,i) * t0( Index(aa,cc,j,k,Nparticles,Nparticles,Nholes)); // Permutation of i,j
                            tau -= v(k,a,c,j) * t0( Index(bb,cc,i,k,Nparticles,Nparticles,Nholes)); // Permutation of a,b
                            tau += v(k,a,c,i) * t0( Index(bb,cc,j,k,Nparticles,Nparticles,Nholes)); // Permutation of a,b,i,j
                        }
                    }

                    // Adding quadratic terms
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){
                            for (int cc=0; cc<Nparticles; cc++){
                                for (int dd=0; dd<Nparticles; dd++){
                                    int c = cc + Nholes; int d = dd + Nholes;

                                    double Qa = 0; double Qb = 0; double Qc = 0; double Qd = 0;

                                    Qa = 0.25*v(k,l,c,d) * t0( Index(cc,dd,i,j,Nparticles,Nparticles,Nholes) ) * t0( Index(aa,bb,k,l,Nparticles,Nparticles,Nholes));

                                    Qb += t0( Index(aa,cc,i,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,bb,l,j,Nparticles,Nparticles,Nholes)); // No permutation
                                    Qb -= t0( Index(aa,cc,j,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,bb,l,i,Nparticles,Nparticles,Nholes)); // Permutation of i,j
                                    Qb -= t0( Index(bb,cc,i,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,aa,l,j,Nparticles,Nparticles,Nholes)); // Permutation of a,b
                                    Qb += t0( Index(bb,cc,j,k,Nparticles,Nparticles,Nholes)) * t0( Index(dd,aa,l,i,Nparticles,Nparticles,Nholes)); // Permutation of a,b,i,j
                                    Qb *= 0.5 * v(k,l,c,d);

                                    Qc -= t0( Index(aa,bb,i,k,Nparticles,Nparticles,Nholes) ) * t0( Index(cc,dd,j,l,Nparticles,Nparticles,Nholes)); // No permutation
                                    Qc += t0( Index(aa,bb,j,k,Nparticles,Nparticles,Nholes) ) * t0( Index(cc,dd,i,l,Nparticles,Nparticles,Nholes)); // Permutation of i,j
                                    Qc *= 0.5 * v(k,l,c,d);

                                    Qd -= t0( Index(bb,dd,k,l,Nparticles,Nparticles,Nholes) ) * t0( Index(aa,cc,i,j,Nparticles,Nparticles,Nholes)); // No permutation
                                    Qd += t0( Index(aa,dd,k,l,Nparticles,Nparticles,Nholes) ) * t0( Index(bb,cc,i,j,Nparticles,Nparticles,Nholes)); // Permutation of a,b
                                    Qd *= 0.5 * v(k,l,c,d);

                                    tau += Qa + Qb + Qc + Qd;
                                }
                            }
                        }
                    }
                    tau = v(a,b,i,j) + tau; //Weighting the iterative scheme

                    t( Index(aa,bb,i,j,Nparticles,Nparticles,Nholes) ) = tau / epsilon(i,j,a,b);
                }
            }
        }
    }
    // Adding weight factor
    if (weight != 0) t = weight*t + (1-weight)*t0;
}

double CCDNaive::epsilon(int i, int j, int a, int b){

    if (BasisNumber == 1) return pabasis.epsilon(i,j,a,b);
    else if (BasisNumber == 2) return elbasis.epsilon(i,j,a,b);
    else {cout << "basis is not defined properly in ccd naive. Epsilon not computed properly" << endl; return 0;}
}

double CCDNaive::v(int p, int q, int r, int s){

    if (BasisNumber == 1) return pabasis.TwoBodyOperator(p,q,r,s);
    else if (BasisNumber == 2) return elbasis.TwoBodyOperator(p,q,r,s);
    else {cout << "basis is not defined properly in ccd naive. TwoBodyOperator not computed properly" << endl; return 0;}
}


double CCDNaive::AbsoluteDifference(double a, double b){
    return sqrt( pow(a-b,2) );
}

int CCDNaive::Index(int p, int q, int r, int s, int Np, int Nq, int Nr){
    // p, q, r, s are the indice and not the state number. i.e. by formalism in this program: aa, bb and not a, b
    // Np, Nq, Nr are the number of indices for each state

    return p + q*Np + r*Np*Nq + s*Np*Nq*Nr;
}
