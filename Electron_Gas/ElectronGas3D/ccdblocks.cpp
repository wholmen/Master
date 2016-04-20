#include "ccdblocks.h"

CCDBlocks::CCDBlocks(){}

CCDBlocks::CCDBlocks(basis_set BASIS){
    basis = BASIS;
    Nholes = basis.Nparticles;
    Nstates = basis.nstates;
    Nparticles = Nstates - Nholes;

    t0 = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);
    t = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);

    UpdateVabcd();
    UpdateVklij();
    UpdateVabij();
    UpdateVklcd();

    I1 = zeros<mat>(Nholes*Nholes, Nholes*Nholes);
    I2 = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);
    I3 = zeros<mat>(Nholes, Nholes);
    I4 = zeros<mat>(Nparticles, Nparticles);
}

double CCDBlocks::Abs(double a, double b){
    return sqrt( pow(a-b, 2) );
}


double CCDBlocks::CCD(){
    int n = 0;
    double error = 0.00001;

    update_t();

    double Eold = 0.0;
    double Enew = CorrolationEnergy();

    while ( n < 10 && Abs(Enew,Eold) > error){

        update_t();
        Eold = Enew;
        Enew = CorrolationEnergy();

        n ++;
    }

    cout << "Coupled Cluster method used " << n << " iterations." << endl;
    return Enew;
}



double CCDBlocks::CorrolationEnergy(){
    double E = 0.0;

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    E += basis.TwoBodyOperator(i,j,a,b) * t(aa+bb*Nparticles, i+j*Nholes);
                }
            }
        }
    }
    return E / 4.0;
}


void CCDBlocks::update_t(){
    t0 = t;
    mat La, Lb, Qa;

    // First calculate La term
    La = Vabcd*t0 / 2.0;


    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Doing a mapping to align diagram Lb
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mat Vtilde = zeros<mat>(Nholes*Nholes, Nholes*Nholes);
    for (int k=0; k<Nholes; k++){
        for (int l=0; l<Nholes; l++){
            for (int i=0; i<Nholes; i++){
                for (int j=0; j<Nholes; j++){

                    Vtilde(k+l*Nholes, i+j*Nholes) = Vklij(i+j*Nholes, k+l*Nholes);
                }
            }
        }
    }
    mat Ttilde = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);

    for (int aa=0; aa<Nparticles; aa++){
        for (int bb=0; bb<Nparticles; bb++){
            for (int k=0; k<Nholes; k++){
                for (int l=0; l<Nholes; l++){

                    Ttilde(k+l*Nholes, aa+bb*Nparticles) = t0(aa+bb*Nparticles, k+l*Nholes);
                }
            }
        }
    }
    mat Lbtilde = 0.5*Vtilde*Ttilde;
    Lb = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);

    for (int aa=0; aa<Nparticles; aa++){
        for (int bb=0; bb<Nparticles; bb++){
            for (int i=0; i<Nholes; i++){
                for (int j=0; j<Nholes; j++){
                    Lb(aa+bb*Nparticles, i+j*Nholes) = Lbtilde(i+j*Nholes, aa+bb*Nparticles);
                }
            }
        }
    }
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Re-mapping finished

    // Following statement will work if mapping was unecessary.
    // Lb = 0.5*Vklij*t0;
    // Statement didn't work. Mapping scheme needed.
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Testing wether a simple transpose of matrix can work
    //Lb = 0.5* Vklij.t() * t0.t();
    //Lb = Lb.t();
    // Working
    // Transposing may be more efficient, but not a general function for all diagrams
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%



    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&
    // Adding the term Qa by a mapping scheme
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&

    Qa = 0.25 * t0 * Vklcd * t0;








    /*
     *This is a start of including Qd
     *
    mat D3tilde = zeros<mat>(Nholes*Nparticles*Nparticles, Nholes);
    mat t1tilde = zeros<mat>(Nholes*Nparticles*Nparticles, Nholes);
    mat Vdtilde = zeros<mat>(Nholes, Nholes*Nparticles*Nparticles);
    mat t2tilde = zeros<mat>(Nholes*Nparticles*Nparticles, Nholes);

    for (int )
    */


    //t = Vabij + 0.5*(La + Lb + Qa);
    UpdateI1();
    mat LI1 = 0.5*t0*I1;

    t = Vabij + 0.5*(La + LI1);

    for (int i=0; i<Nholes; i++){
        for (int j=0; j<Nholes; j++){
            for (int aa=0; aa<Nparticles; aa++){
                for (int bb=0; bb<Nparticles; bb++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    t(aa+bb*Nparticles, i+j*Nholes) /= basis.epsilon(i,j,a,b);
                }
            }
        }
    }
}



void CCDBlocks::UpdateVabcd(){
    Vabcd = zeros<mat>(Nparticles*Nparticles, Nparticles*Nparticles); // Interaction matrix based on La term

    for (int aa=0; aa<Nparticles; aa++){
        for (int bb=0; bb<Nparticles; bb++){
            for (int cc=0; cc<Nparticles; cc++){
                for (int dd=0; dd<Nparticles; dd++){
                    int a = aa + Nholes; int b = bb + Nholes;
                    int c = cc + Nholes; int d = dd + Nholes;

                    Vabcd(aa + bb*Nparticles, cc + dd*Nparticles) = basis.TwoBodyOperator(a,b,c,d);
                }
            }
        }
    }
}

void CCDBlocks::UpdateVabij(){
    Vabij = zeros<mat>(Nparticles*Nparticles, Nholes*Nholes);

    for (int aa=0; aa<Nparticles; aa++){
        for (int bb=0; bb<Nparticles; bb++){
            for (int i=0; i<Nholes; i++){
                for (int j=0; j<Nholes; j++){
                    int a = aa + Nholes; int b = bb + Nholes;

                    Vabij(aa+bb*Nparticles, i+j*Nholes) = basis.TwoBodyOperator(a,b,i,j);
                }
            }
        }
    }
}

void CCDBlocks::UpdateVklij(){
    Vklij = zeros<mat>(Nholes*Nholes, Nholes*Nholes);

    for (int k=0; k<Nholes; k++){
        for (int l=0; l<Nholes; l++){
            for (int i=0; i<Nholes; i++){
                for (int j=0; j<Nholes; j++){

                    Vklij(k+l*Nholes, i+j*Nholes) = basis.TwoBodyOperator(k,l,i,j);
                }
            }
        }
    }
}

void CCDBlocks::UpdateVklcd(){
    Vklcd = zeros<mat>(Nholes*Nholes, Nparticles*Nparticles);

    for (int k=0; k<Nholes; k++){
        for (int l=0; l<Nholes; l++){
            for (int cc=0; cc<Nparticles; cc++){
                for (int dd=0; dd<Nparticles; dd++){
                    int c = cc + Nholes; int d = dd + Nholes;

                    Vklcd(k+l*Nholes, cc+dd*Nparticles) = basis.TwoBodyOperator(k,l,c,d);
                }
            }
        }
    }
}

void CCDBlocks::UpdateI1(){

    I1 = Vklij + 0.5 * Vklcd*t;
}
















