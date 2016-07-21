#include "mbptnaive.h"

MBPTNaive::MBPTNaive()
{

}

MBPTNaive::MBPTNaive(PairingBasis BASIS){
    basis = BASIS;

    // Calculating important variables
    Nholes = basis.Nholes; Nholes2 = Nholes*Nholes; Nholes3 = Nholes2*Nholes;
    Nstates = basis.Nstates;
    Nparticles = Nstates - Nholes; Nparticles2 = Nparticles*Nparticles; Nparticles3 = Nparticles2*Nparticles;

}


double MBPTNaive::MBPT2_MHJ(){

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


double MBPTNaive::MBPT2(){

    double D1 = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){

            for (int a=Nholes; a<Nstates; a++){
                for (int b=Nholes; b<Nstates; b++){

                    D1 += basis.TwoBodyOperator(i,j,a,b)*basis.TwoBodyOperator(a,b,i,j) / basis.epsilon(i,j,a,b);
                }
            }
        }
    }
    D1 /= 4.0;
    return D1;
}

double MBPTNaive::MBPT3(){

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

double MBPTNaive::MBPT4(){

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
                                    D35 += basis.TwoBodyOperator(a,c,i,j)*basis.TwoBodyOperator(i,j,a,b)*basis.TwoBodyOperator(b,d,k,l)*basis.TwoBodyOperator(k,l,c,d) / (basis.epsilon(k,l,c,d)*basis.epsilon(i,j,a,b)*basis.epsilon4(i,j,k,l,a,b,c,d));
                                    D36 += basis.TwoBodyOperator(i,j,a,b)*basis.TwoBodyOperator(k,l,c,d)*basis.TwoBodyOperator(a,b,k,l)*basis.TwoBodyOperator(c,d,i,j) / (basis.epsilon(i,j,a,b)*basis.epsilon(i,j,c,d)*basis.epsilon4(i,j,k,l,a,b,c,d));
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
