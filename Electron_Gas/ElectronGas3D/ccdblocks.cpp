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
    UpdateVkbcj();

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

    bool Intermediates = true;



    if (Intermediates){

        t = Vabij + 0.5*(DiagramLa() + DiagramI1() + DiagramI2() );
    }
    else{

        t = Vabij + 0.5*(DiagramLa() + DiagramLb() + DiagramQa());
    }


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

mat CCDBlocks::DiagramLa(){
    // First calculate La term
    mat La;
    La = Vabcd*t0 / 2.0;
    return La;
}

mat CCDBlocks::DiagramLb(){
    // Function calculates diagram Lb

    mat Vtilde = Realign(Vklij, 2,Nholes, 3,Nholes, 0,Nholes, 1,Nholes);

    mat Ttilde = Realign(t0, 2,Nparticles, 3,Nparticles, 0,Nholes, 1,Nholes);

    mat Lbtilde = 0.5*Vtilde*Ttilde;

    mat Lb = Realign(Lbtilde, 2,Nholes, 3,Nholes, 0, Nparticles, 1, Nparticles);

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Testing wether a simple transpose of matrix can work
    //Lb = 0.5* Vklij.t() * t0.t();
    //Lb = Lb.t();
    // Working
    // Transposing may be more efficient, but not a general function for all diagrams
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return Lb;
}

mat CCDBlocks::DiagramQa(){
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&
    // Adding the term Qa by a mapping scheme
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&

    mat Qa = 0.25 * t0 * Vklcd * t0;
    return Qa;
}


mat CCDBlocks::DiagramI1(){
    UpdateI1();
    return 0.5*t0*I1;
}

mat CCDBlocks::DiagramI2(){
    UpdateI2();

    mat Ttilde = Realign( t0, 0,Nparticles, 2,Nparticles, 1,Nholes, 3,Nholes);

    mat Itilde = Realign( I2, 1,Nholes, 2,Nparticles, 0,Nparticles, 3,Nholes);

    return Realign( Ttilde*Itilde, 0,Nparticles, 2,Nholes, 1,Nparticles, 3,Nholes);
}






// Functions to set up two body interaction matrices and intermediates

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

void CCDBlocks::UpdateVkbcj(){
    Vkbcj = zeros<mat>(Nholes*Nparticles, Nholes*Nparticles);

    for (int k=0; k<Nholes; k++){
        for (int bb=0; bb<Nparticles; bb++){
            for (int cc=0; cc<Nparticles; cc++){
                for (int j=0; j<Nholes; j++){
                    int b = bb + Nholes; int c = cc + Nholes;

                    Vkbcj(bb+k*Nparticles, cc+j*Nparticles) = basis.TwoBodyOperator(k,b,c,j);
                }
            }
        }
    }
}

void CCDBlocks::UpdateI1(){
    I1 = Vklij + 0.5 * Vklcd*t;
}

void CCDBlocks::UpdateI2(){
    // I2 consist of 4 parts because of

    mat Vtilde = Realign( Vklcd, 0,Nholes, 2,Nholes, 1,Nparticles, 3,Nparticles);

    mat Ttilde = Realign( t0, 1,Nparticles, 2,Nparticles, 0,Nholes, 3,Nholes);

    I2 = Vkbcj + 0.5 * Realign( Vtilde*Ttilde , 0,Nholes, 2,Nparticles, 1,Nparticles, 3,Nholes);
}





mat CCDBlocks::Realign(mat A, int a1, int Na, int b1, int Nb, int c1, int Nc, int d1, int Nd){
    // A function set to realign a matrix A.
    // This function will restructure the size of matrix A

    // One uses the function by sending in the new placements of indices
    // Example
    //         a1 = 2, b1 = 1, c1 = 3, d1 = 0
    //         A (a+b*Na, c+d*Nc)  ->  A (d+b*Nd, a+c*Na)

    vec e = zeros<vec>(4);
    vec eN = zeros<vec>(4);

    eN(a1) = Na; eN(b1) = Nb; eN(c1) = Nc; eN(d1) = Nd;

    mat Anew = zeros<mat>( eN(0)*eN(1), eN(2)*eN(3) );

    for (int a=0; a<Na; a++){
        for (int b=0; b<Nb; b++){
            for (int c=0; c<Nc; c++){
                for (int d=0; d<Nd; d++){

                    e(a1)  =  a; e(b1) = b; e(c1) = c; e(d1) = d;
                    eN(a1) = Na; eN(b1) = Nb; eN(c1) = Nc; eN(d1) = Nd;

                    Anew( e(0)+e(1)*eN(0), e(2)+e(3)*eN(2) ) = A ( a+b*Na, c+d*Nc );
                }
            }
        }
    }
    return Anew;
}

















