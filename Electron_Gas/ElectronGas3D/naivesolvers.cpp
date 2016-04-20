#include "naivesolvers.h"

NaiveSolvers::NaiveSolvers()
{

}

NaiveSolvers::NaiveSolvers(basis_set BASIS){
    basis = BASIS;
    Nholes = basis.Nparticles;
    Nstates = basis.nstates;
    Nparticles = Nstates - Nholes;
}


double NaiveSolvers::CCD(){
    double Ecc_old = 0.0; double Ecc_new = 0.0;
    int n = 0;

    update_t0();
    Ecc_old = dE_CCD();

    update_t(2);
    Ecc_new = dE_CCD();

    while ( n < 10 && sqrt( pow(Ecc_new - Ecc_old,2) ) > 0.00001){

        Ecc_old = Ecc_new;
        update_t(2);
        Ecc_new = dE_CCD();
        n ++;
    }
    cout << "Coupled Cluster method used " << n << " iterations." << endl;
    return Ecc_new;
}


double NaiveSolvers::dE_CCD(){
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

void NaiveSolvers::update_t0(){
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


void NaiveSolvers::update_t(int degree){
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

                            tau += 0.5 * basis.TwoBodyOperator(a,b,c,d) * t0(i+j*Nholes, cc+dd*Nparticles);
                        }
                    }

                    // Adding Lb
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            tau += 0.5 * basis.TwoBodyOperator(k,l,i,j) * t0(k+l*Nholes, aa+bb*Nparticles);
                        }
                    }
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
                    if (degree == 2){
                        for (int k=0; k<Nholes; k++){
                            for (int l=0; l<Nholes; l++){
                                for (int cc=0; cc<Nparticles; cc++){
                                    for (int dd=0; dd<Nparticles; dd++){
                                        int c = cc + Nholes; int d = dd + Nholes;

                                        double Qa = 0; double Qb = 0; double Qc = 0; double Qd = 0;

                                        Qa = 0.25*basis.TwoBodyOperator(k,l,c,d) * t0(i+j*Nholes, cc+dd*Nparticles) * t0(k+l*Nholes, aa+bb*Nparticles);
/*
                                        Qb += t0(i+k*Nholes, aa+cc*Nparticles) * t0(l+j*Nholes, dd+bb*Nparticles); // No permutation
                                        Qb -= t0(j+k*Nholes, aa+cc*Nparticles) * t0(l+i*Nholes, dd+bb*Nparticles); // Permutation of i,j
                                        Qb -= t0(i+k*Nholes, bb+cc*Nparticles) * t0(l+j*Nholes, dd+aa*Nparticles); // Permutation of a,b
                                        Qb += t0(j+k*Nholes, bb+cc*Nparticles) * t0(l+i*Nholes, dd+aa*Nparticles); // Permutation of a,b,i,j
                                        Qb *= 0.5 * basis.TwoBodyOperator(k,l,c,d);

                                        Qc -= t0(i+k*Nholes, aa+bb*Nparticles) * t0(j+l*Nholes, cc+dd*Nparticles); // No permutation
                                        Qc += t0(j+k*Nholes, aa+bb*Nparticles) * t0(i+l*Nholes, cc+dd*Nparticles); // Permutation of i,j
                                        Qc *= 0.5 * basis.TwoBodyOperator(k,l,c,d);

                                        Qd -= t0(k+l*Nholes, bb+dd*Nparticles) * t0(i+j*Nholes, aa+cc*Nparticles); // No permutation
                                        Qd += t0(k+l*Nholes, aa+dd*Nparticles) * t0(i+j*Nholes, bb+cc*Nparticles); // Permutation of a,b
                                        Qd *= 0.5 * basis.TwoBodyOperator(k,l,c,d);
*/
                                        tau += Qa + Qb + Qc + Qd;
                                    }
                                }
                            }
                        }
                    }

                    tau = basis.TwoBodyOperator(a,b,i,j) + 0.5*tau; //Weighting the iterative scheme

                    t(i + j*Nholes, aa + bb*Nparticles) = tau / basis.epsilon(i,j,a,b);
                }
            }
        }
    }

}


double NaiveSolvers::MBPT2_MHJ(){

    mat vhhpp = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);
    mat vpphh = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            for (int a=Nholes; a<Nstates; a++){
                for (int b=Nholes; b<Nstates; b++){

                    int aa = a - Nholes; int bb = b - Nholes;

                    vhhpp(i + j*Nholes, aa + bb*Nparticles) = basis.TwoBodyOperator(i,j,a,b);
                    vpphh(aa + bb*Nparticles, i + j*Nholes) = basis.TwoBodyOperator(a,b,i,j) / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
    return 0.25 * trace(vhhpp*vpphh);

}


double NaiveSolvers::MBPT2(){

    double D1 = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            for (int a=Nholes; a<Nstates; a++){
                for (int b=Nholes; b<Nstates; b++){
                    //cout << "i: " << i << "  , j: " << j << "  , a: " << a << "  , b: " << b << endl;
                    D1 += basis.TwoBodyOperator(i,j,a,b)*basis.TwoBodyOperator(a,b,i,j) / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
    D1 /= 4;
    return D1;
}

double NaiveSolvers::MBPT3(){

    double D1 = MBPT2();

    double D4 = 0; double D5 = 0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            for (int a=Nholes; a<Nstates; a++){
                for (int b=Nholes; b<Nstates; b++){

                    for (int c=Nholes; c<Nstates; c++){
                        for (int d=Nholes; d<Nstates; d++){

                            D4 += basis.TwoBodyOperator(c,d,i,j)*basis.TwoBodyOperator(a,b,c,d)*basis.TwoBodyOperator(i,j,a,b) / basis.epsilon(i,j,c,d) / basis.epsilon(i,j,a,b);
                        }
                    }
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){

                            D5 += basis.TwoBodyOperator(a,b,k,l)*basis.TwoBodyOperator(k,l,i,j)*basis.TwoBodyOperator(i,j,a,b) / basis.epsilon(k,l,a,b) / basis.epsilon(i,j,a,b);
                        }
                    }
                }
            }
        }
    }
    D4 /= 8.0;
    D5 /= 8.0;
    return D1 + D4 + D5;
}

double NaiveSolvers::MBPT4(){

    double mbpt3 = MBPT3();
    double D5, D6, D14, D15, D34, D35, D36, D37, D38, D39, D40;
    D5 = D6 = D14 = D15 = D34 = D35 = D36 = D37 = D38 = D39 = D40 = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            for (int a=Nholes; a<Nstates; a++){
                for (int b=Nholes; b<Nstates; b++){

                    for (int c=Nholes; c<Nstates; c++){
                        for (int d=Nholes; d<Nstates; d++){
                            for (int k=0; k<Nholes; k++){
                                for (int l=0; l<Nholes; l++){

                                    D5  += basis.TwoBodyOperator(c,d,k,l)*basis.TwoBodyOperator(k,l,i,j)*basis.TwoBodyOperator(a,b,c,d)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(k,l,c,d)*basis.epsilon(i,j,c,d)*basis.epsilon(i,j,a,b));
                                    D6  += basis.TwoBodyOperator(c,d,k,l)*basis.TwoBodyOperator(k,l,i,j)*basis.TwoBodyOperator(a,b,c,d)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(k,l,c,d)*basis.epsilon(k,l,a,b)*basis.epsilon(i,j,a,b));

                                    D34 += basis.TwoBodyOperator(a,b,i,k)*basis.TwoBodyOperator(i,j,a,b)*basis.TwoBodyOperator(c,d,j,l)*basis.TwoBodyOperator(k,l,c,d) / (basis.epsilon(j,l,c,d)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D35 += basis.TwoBodyOperator(a,c,i,j)*basis.TwoBodyOperator(i,j,a,b)*basis.TwoBodyOperator(c,d,k,l)*basis.TwoBodyOperator(k,l,c,d) / (basis.epsilon(k,l,c,d)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D36 += basis.TwoBodyOperator(b,d,i,l)*basis.TwoBodyOperator(a,c,j,k)*basis.TwoBodyOperator(j,k,b,d)*basis.TwoBodyOperator(a,c,i,l) / (basis.epsilon(i,l,b,d)*basis.epsilon(i,l,a,c)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D37 += basis.TwoBodyOperator(a,b,k,l)*basis.TwoBodyOperator(c,d,i,j)*basis.TwoBodyOperator(k,l,c,d)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(k,l,a,b)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D38 += basis.TwoBodyOperator(a,c,i,k)*basis.TwoBodyOperator(k,l,c,d)*basis.TwoBodyOperator(d,b,l,j)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(k,i,a,c)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D39 += basis.TwoBodyOperator(a,c,i,j)*basis.TwoBodyOperator(b,d,k,l)*basis.TwoBodyOperator(k,l,c,d)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(i,j,a,c)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D40 += basis.TwoBodyOperator(a,b,i,k)*basis.TwoBodyOperator(c,d,j,l)*basis.TwoBodyOperator(k,l,c,d)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(i,k,a,b)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                }
                            }
                        }
                    }
                    for (int c=Nholes; c<Nstates; c++){
                        for (int d=Nholes; d<Nstates; d++){
                            for (int e=Nholes; e<Nstates; e++){
                                for (int f=Nholes; f<Nstates; f++){

                                    D14 += basis.TwoBodyOperator(e,f,i,j)*basis.TwoBodyOperator(c,d,e,f)*basis.TwoBodyOperator(a,b,c,d)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(i,j,e,f)*basis.epsilon(i,j,c,d)*basis.epsilon(i,j,a,b));
                                }
                            }
                        }
                    }
                    for (int k=0; k<Nholes; k++){
                        for (int l=0; l<Nholes; l++){
                            for (int m=0; m<Nholes; m++){
                                for (int n=0; n<Nholes; n++){
                                    D15 += basis.TwoBodyOperator(a,b,m,n)*basis.TwoBodyOperator(m,n,k,l)*basis.TwoBodyOperator(k,l,i,j)*basis.TwoBodyOperator(i,j,a,b) / (basis.epsilon(m,n,a,b)*basis.epsilon(k,l,a,b)*basis.epsilon(i,j,a,b));

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    D5 /= 16.0; D6 /= 16.0; D14 /= 16.0; D15 /= 16.0;
    D34 /= -4.0; D35 /= -4.0; D36 /= 16.0; D37 /= 16.0; D38 /= 1.0; D39 /= -4.0; D40 /= -4.0;

    return mbpt3 + D5 + D6 + D14 + D15 + D34 + D35 + D36 + D37 + D38 + D39 + D40;
}



