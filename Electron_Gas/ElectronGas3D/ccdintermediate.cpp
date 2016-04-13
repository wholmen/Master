#include "ccdintermediate.h"

CCDIntermediate::CCDIntermediate(){}

CCDIntermediate::CCDIntermediate(basis_set BASIS){
    basis = BASIS;
    Nholes = basis.Nparticles;
    Nstates = basis.nstates;
    Nparticles = Nstates - Nholes;

    I1 = zeros<mat>(Nholes*Nholes, Nholes*Nholes);
    I2 = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);
    I3 = zeros<mat>(Nholes, Nholes);
    I4 = zeros<mat>(Nparticles, Nparticles);
}


double CCDIntermediate::CCD(){
    double Ecc_old = 0.0; double Ecc_new = 0.0;
    int n = 0;

    update_t0();
    Ecc_old = dE_CCD();

    update_t();
    Ecc_new = dE_CCD();

    while ( n < 10 && sqrt( pow(Ecc_new - Ecc_old,2) ) > 0.00001){

        Ecc_old = Ecc_new;
        update_t();
        Ecc_new = dE_CCD();
        n ++;
    }
    cout << "Coupled Cluster method used " << n << " iterations." << endl;
    return Ecc_new;
}



double CCDIntermediate::dE_CCD(){
    double E = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += basis.TwoBodyOperator(i,j,a,b) * t(i + j*Nholes, aa + bb*Nparticles);
                }
            }
        }
    }
    return E * 0.25;
}

void CCDIntermediate::update_t0(){
    // Setting up the first amplitudes based on a guess based on mbpt
    t0 = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    t0(i + j*Nholes, aa + bb*Nparticles) = basis.TwoBodyOperator(a,b,i,j) / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
    t = t0;
}

void CCDIntermediate::update_t(){
    t0 = t;

    UpdateI1();
    UpdateI2();
    UpdateI3();
    UpdateI4();

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
                            term += basis.TwoBodyOperator(a,b,c,d) * t0(i+j*Nholes, cc+dd*Nparticles);
                        }
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            term += I1(i+j*Nholes, k+l*Nholes) * t0(k+l*Nholes, aa+bb*Nparticles);
                        }
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        for (int cc=0; cc<Nparticles; cc++){

                            term += I2(j+k*Nholes, bb+cc*Nparticles) * t0(i+k*Nholes, aa+cc*Nparticles); // No permutation
                            term -= I2(i+k*Nholes, bb+cc*Nparticles) * t0(j+k*Nholes, aa+cc*Nparticles); // Permutation i,j
                            term -= I2(j+k*Nholes, aa+cc*Nparticles) * t0(i+k*Nholes, bb+cc*Nparticles); // Permutation a,b
                            term += I2(i+k*Nholes, aa+cc*Nparticles) * t0(j+k*Nholes, bb+cc*Nparticles); // Permutation a,b,i,j
                        }
                    }
                    tau += term;

                    term = 0;
                    for (int k=0; k<Nholes; k++){
                        term -= I3(j,k) * t0(i+k*Nholes, aa+bb*Nparticles); // No permutation
                        term += I3(i,k) * t0(j+k*Nholes, aa+bb*Nparticles); // Permutation i,j
                    }
                    tau += 0.5*term;

                    term = 0;
                    for (int cc=0; cc<Nparticles; cc++){
                        term -= I4(bb,cc) * t0(i+j*Nholes, aa+cc*Nparticles); // No permutation
                        term += I4(aa,cc) * t0(i+j*Nholes, bb+cc*Nparticles); // Permutation a,b
                    }
                    tau += 0.5*term;

                    tau = basis.TwoBodyOperator(a,b,i,j) + 0.5*tau; // Weighting the iterative scheme

                    t(i+j*Nholes, aa+bb*Nparticles) = tau / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
}




void CCDIntermediate::UpdateI1(){
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
                            I1(i+j*Nholes, k+l*Nholes) += basis.TwoBodyOperator(k,l,c,d) * t(i+j*Nholes, cc+dd*Nparticles);
                        }
                    }
                    I1(i+j*Nholes, k+l*Nholes) = basis.TwoBodyOperator(k,l,i,j) + 0.5*I1(i+j*Nholes, k+l*Nholes);
                }
            }
        }
    }
}

void CCDIntermediate::UpdateI2(){
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
                            I2(j+k*Nholes, bb+cc*Nparticles) += basis.TwoBodyOperator(k,l,c,d) * t(l+j*Nholes, dd+bb*Nparticles);
                        }
                    }
                    I2(j+k*Nholes, bb+cc*Nparticles) = basis.TwoBodyOperator(k,b,c,j) + 0.5*I2(j+k*Nholes, bb+cc*Nparticles);
                }
            }
        }
    }
}

void CCDIntermediate::UpdateI3(){
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

                        I3(j,k) += basis.TwoBodyOperator(k,l,c,d) * t(j+l*Nholes, cc+dd*Nparticles);
                    }
                }
            }
        }
    }
}

void CCDIntermediate::UpdateI4(){
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

                        I4(bb,cc) += basis.TwoBodyOperator(k,l,c,d) * t(k+l*Nholes, bb+dd*Nparticles);
                    }
                }
            }
        }
    }
}

























